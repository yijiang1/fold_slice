/*!
 * \file
 * Read Object Data
 *
 * This file contains the main mex code for reading object data into MATLAB memory.
 */

#include <string>
#include <memory>
#include <chrono>
#include <cassert>
#include <array>
#include <vector>
#include <cstring>
#include <cerrno>
#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <unistd.h>
#include <functional>
#include <algorithm>
#include <sys/mman.h>
#include <sys/wait.h>
#include "precision.h"
#include "hdf5_helper.h"
#include "mex.h"
#include "read_object_data.h"
#include "debug_helper.h"
#include "multi_processing.h"

namespace {

    constexpr int min_dims = 2;                 //!< Minimum number of dimensions (row, col)
    constexpr int max_dims = 4;                 //!< Maximum number of dimensions (nslices, nmodes, row, col)

    #include "mex_helper.h"

    /*!
     * \brief Read error exception class
     */
    struct read_error : public std::exception {
        bool converted = false;    //!< has the error message been converted (to include 1-based file index)
        uint32_t file_index;       //!< file index (starting from 0)
        std::string error_message; //!< error message (converted if converted is true)

        /*!
         * \brief Constructor
         * \param index 0-based file index
         * \param message unconverted error message
         */
        read_error(uint32_t index, const std::string &message)
            : file_index(index), error_message(message)
        {}

        /*!
         * \brief Retrieve converted error message
         * This method first attempts to convert the error message by appending the 1-based file index (MATLAB indexing convention).
         * If this fails, the unconverted error message is returned.
         * \return converted error message
         */
        const char* what() const noexcept override
        {
            try {
                if (! converted) {
                    std::ostringstream oss;
                    oss << error_message << " (file " << (file_index+1) << ")";
                    const_cast<std::string&>(error_message) = oss.str();
                    const_cast<bool&>(converted) = true;
                }
            } catch (...) {}
            return error_message.c_str();
        }
    };

    /*!
     * \brief Error collector class
     * The class either keeps a vector to collect file indices for files that could not be read, or it throws an exception, depending on the in_collect_mode flag.
     * The class provides a pointer to the singleton instance.
     */
    struct error_collector {
        static error_collector *singleton;  //!< static pointer to make the singleton instance available
        bool in_collect_mode = false;       //!< is the error method collecting or throwing errors?
        std::vector<uint32_t> bad_index;    //!< vector of bad file indices

        /*!
         * \brief Constructor
         * \param collect_mode true if errors should be collected instead of thrown
         */
        explicit error_collector(bool collect_mode)
            : in_collect_mode(collect_mode)
        {
            if (! singleton)
                singleton = this;
            else
                throw std::runtime_error("internal error - multiple error_collector instances!");
        }

        /*!
         * \brief Destructor
         * Zero out the singleton pointer
         */
        ~error_collector() noexcept
        {
            singleton = nullptr;
        }

        /*!
         * \brief Get pointer to singleton instance
         */
        static error_collector& get() noexcept
        {
            return *singleton;
        }

        /*!
         * \brief Collect or throw error
         * \param ex read error exception that would be thrown if not in collector mode
         */
        void error(read_error &&ex)
        {
            if (! in_collect_mode)
                throw ex;
            bad_index.push_back(ex.file_index);
            DEBUG {
                OUT << "error reading file " << ex.file_index << ": " << ex.error_message << std::endl;
            }
        }

        /*!
         * \brief Add index to bad indices vector
         * \param index index of file that can not be read
         */
        void add(uint32_t index) noexcept
        {
            bad_index.push_back(index);
        }
    };

    error_collector* error_collector::singleton = nullptr;  //!< static error_collectr singleton pointer

    /*!
     * \brief Adapt data space layout
     *
     * Adapt the memory and file dataset layouts by selecting apropriate hyperslabs
     *
     * \param mspace memory data space layout
     * \param fspace file dataset layout
     * \param mdim memory space dimensions
     * \param fdim file dataset dimensions
     * \param doffset Start of the dimsensions present in the file data
     */
    void adapt_spaces(hdf5::dataspace &mspace, hdf5::dataspace &fspace, const hsize_t mdim[max_dims], const hsize_t fdim[max_dims], int doffset)
    {
        hsize_t mstart[max_dims] = {0};
        hsize_t fstart[max_dims] = {0};
        hsize_t mcount[max_dims] = { mdim[0], mdim[1], mdim[2], mdim[3] };
        hsize_t fcount[max_dims] = { fdim[0], fdim[1], fdim[2], fdim[3] };
        bool mselect = false;
        bool fselect = false;

        for (int i=2; i<=3; i++) {
            if (mdim[i] < fdim[i]) {
                fstart[i] = (fdim[i] - mdim[i]) / 2;
                fcount[i] = mdim[i];
                fselect = true;
            } else if (mdim[i] > fdim[i]) {
                mstart[i] = (mdim[i] - fdim[i]) / 2;
                mcount[i] = fdim[i];
                mselect = true;
            }
        }
        if (mselect) {
            DEBUG {
                OUT << "memory subspace start=" << mstart[0] << ',' << mstart[1] << ',' << mstart[2] << ',' << mstart[3] << '/' << doffset << " count=" << mcount[0] << ',' << mcount[1] << ',' << mcount[2] << ',' << mcount[3] << '/' << doffset << "\nfile space is " << fdim[0] << ',' << fdim[1] << ',' << fdim[2] << ',' << fdim[3] << std::endl;
            }
            if (H5Sselect_hyperslab(mspace, H5S_SELECT_SET, &mstart[doffset], nullptr, &mcount[doffset], nullptr) < 0)
                throw hdf5::exception("unable to set hdf5 memory space");
        }
        if (fselect) {
            DEBUG {
                OUT << "file subspace start=" << fstart[0] << ',' << fstart[1] << ',' << fstart[2] << ',' << fstart[3] << '/' << doffset << " count=" << fcount[0] << ',' << fcount[1] << ',' << fcount[2] << ',' << fcount[3] << '/' << doffset << "\nmem space is " << mdim[0] << ',' << mdim[1] << ',' << mdim[2] << ',' << mdim[3] << std::endl;
            }
            if (H5Sselect_hyperslab(fspace, H5S_SELECT_SET, &fstart[doffset], nullptr, &fcount[doffset], nullptr) < 0)
                throw hdf5::exception("unable to set hdf5 dataset space");
        }
    }

    /*!
     * \brief Read object from file into buffer
     * \param file HDF5 file
     * \param object_path path to object dataset
     * \param mdims result matrix dimensions (with size 5)
     * \param buf pointer to parent acessible buffer
     * \tparam f_type buffer element type
     */
    template <typename f_type>
    void read_complex_matrix(hdf5::file &file, const std::string &object_path,
                             const std::vector<mwSize> &mdims, typename mx_trait<f_type>::complex_type *buf)
    {
        hsize_t fdims[max_dims] = { 1, 1, 0, 0 }; // (nslices, nmodes, nrows, ncols)
        hsize_t bdims[max_dims] = { 1, 1, 0, 0 }; // file and buffer dimensions
        unsigned long nelements = 0;
        int dim_offset = 0;
        int ndims = 0;
        int mdim_offset = 0;
        for (; (mdim_offset < max_dims-2) && (mdims[mdim_offset] == 1); mdim_offset++);

        // Create complex type
        hdf5::h5type complex_type(hdf5::create_complex_type<f_type>());
        // Open dataset
        hdf5::dataset dataset(H5Dopen(file, object_path.c_str(), H5P_DEFAULT));
        if (! dataset.valid())
            throw hdf5::exception("unable to open dataset");
        hdf5::dataspace dataspace(H5Dget_space(dataset));
        if (! dataspace.valid())
            throw hdf5::exception("unable to open data space for dataset");
        ndims = H5Sget_simple_extent_ndims(dataspace);
        if (ndims < min_dims || ndims > max_dims)
            throw hdf5::exception("wrong number of dimensions in dataset");
        dim_offset = max_dims - ndims;
        if (H5Sget_simple_extent_dims(dataspace, &fdims[dim_offset], NULL) != ndims)
            throw hdf5::exception("dimension error");
        for (int i=max_dims-2; i; i--) {
            if (fdims[i-1] != mdims[i]) {
                if (i != mdim_offset)
                    throw hdf5::exception("dimension mismatch");
            }
        }
        // Create memory space and adapt data space
        bdims[0] = fdims[0]; // nslices
        bdims[1] = fdims[1]; // nmodes
        bdims[2] = mdims[3]; // nrows wanted
        bdims[3] = mdims[4]; // ncols wanted
        DEBUG {
            OUT << "bdims(" << bdims[0] << ',' << bdims[1] << ',' << bdims[2] << '/' << fdims[2] << ',' << bdims[3] << '/' << fdims[3] << ')' << dim_offset << std::endl;
        }
        hdf5::dataspace memspace(H5Screate_simple(ndims, &bdims[dim_offset], nullptr));
        if (! memspace.valid())
            throw hdf5::exception("unable to create memory data space");
        if ((bdims[2] != fdims[2]) || (bdims[3] != fdims[3]))
            adapt_spaces(memspace, dataspace, bdims, fdims, dim_offset);
        // Read dataset
        if (H5Dread(dataset, complex_type, memspace, dataspace, H5P_DEFAULT, buf) < 0)
            throw hdf5::exception("unable to read dataset");

        DEBUG {
            int last = bdims[max_dims-2] * bdims[max_dims-1] - 1;
            OUT << "read finished - [0]=(" << buf[0].real << ',' << buf[0].imag << ") [" << last << "]=(" << buf[last].real << ',' << buf[last].imag << ')' << std::endl;
        }
    }

    /*!
     * \brief Open HDF5 file
     * \param file_path path to HDF5 file
     * \return HDF5 file object id
     */
    hid_t open_file(const std::string &file_path)
    {
        hdf5::file file(H5Fopen(file_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
        if (! file.valid())
            throw hdf5::exception("unable to open file");
        return file.grab();
    }

    /*!
     * \brief Nonzero dimension size?
     * \param dims size of twodimensional object
     * \tparam T type of dimension size
     * \return are dimension sizes positive?
     */
    template <typename T>
    bool dims_ok (const T *dims) noexcept
    {
        return dims[0] > 0 || dims[1] > 0;
    }

    /*!
     * \brief Drop leading dimensions of size 1
     * \param dim dimsension sizes
     * \param max maximum dimensions
     * \tparam T type of dimension size
     * \return number of dimsension of size bigger than 1
     */
    template <typename T>
    int collapse_dims (const std::vector<T> &dim, int max) noexcept
    {
        int i=0;
        for (; i<max-2; i++) {
            if (dim[i] != 1)
                break;
        }
        return max-i;
    }

    /*!
     * \brief Get HDF5 object size
     *
     * \param file HDF5 file
     * \param dataset_path object dataset path
     * \param dims filled with object dimensions
     */
    void get_object_size(hdf5::file &file, const std::string &dataset_path, hsize_t dims[max_dims])
    {
        int dim_offset = 0;
        int ndims;

        hdf5::dataset dataset(H5Dopen(file, dataset_path.c_str(), H5P_DEFAULT));
        if (! dataset.valid())
            throw hdf5::exception("unable to open dataset");
        hdf5::dataspace dataspace(H5Dget_space(dataset));
        if (! dataspace.valid())
            throw hdf5::exception("unable to get data space");
        ndims = H5Sget_simple_extent_ndims(dataspace);
        if (ndims < min_dims || ndims > max_dims)
            throw hdf5::exception("wrong number of dimensions in dataset");
        dim_offset = max_dims - ndims;
        if (H5Sget_simple_extent_dims(dataspace, &dims[dim_offset], NULL) != ndims)
            throw hdf5::exception("unable to get dimension sizes");
    }

    /*!
     * \brief Read in objects with multiple processes
     *
     * \param nprocs number of read processes
     * \param dims desired result dimensions (or [0, 0] to get the default behaviour)
     * \param object_path path to hdf5 object within the files
     * \param file_paths paths to hdf5 objects
     * \param r_stat fill in if not null
     * \tparam f_type float or double result type
     * \return MATLAB result array
     */
    template <typename f_type>
    mxArray* read_objects_parallel(int nprocs, std::array<int64_t, 2> &dims,
                                   const std::string &object_path, const std::vector<std::string> &file_paths,
                                   data_prep::read_stat *r_stat)
    {
        try {
            std::vector<mwSize> mdims({ 1, 1, 1, 0, 0 });   // extended matrix dimensions
            hsize_t odims[max_dims] = { 1, 1, 0, 0 };       // extended object dimensions (nobjects, nslices, nmodes, nrows, ncols)

            // Assign matrix dimension values
            {
                hdf5::file file(open_file(file_paths[0]));
                get_object_size(file, object_path, odims);
            }
            for (unsigned int i=1; i<max_dims-1; i++) {
                mdims[i] = (mwSize)odims[i-1];
            }
            if (!dims_ok(&dims[0])) {
                mdims[max_dims-1] = odims[max_dims-2];
                mdims[max_dims] = odims[max_dims-1];
            } else {
                mdims[max_dims-1] = dims[0];
                mdims[max_dims] = dims[1];
            }
            int nmdims = collapse_dims(mdims, max_dims+1);
            if (file_paths.size() > 1) {
                mdims[max_dims - nmdims] = (mwSize)file_paths.size();
                nmdims++;
            }

            std::vector<mwSize> matrix_dims(mdims);         // matrix dimensions
            matrix_dims.erase(matrix_dims.begin(), matrix_dims.begin() + max_dims + 1 - nmdims); // cut away extension dimensions

            // Create space for MATLAB array
            mx_ptr<mxArray> matrix;
            {
                std::vector<mwSize> rdims(matrix_dims);
                std::reverse(rdims.begin(), rdims.end());
                matrix.reset(mxCreateNumericArray(nmdims, &rdims[0], mx_trait<f_type>::class_id, mxCOMPLEX));
            }
            if (! matrix.get())
                throw std::runtime_error("matrix creation failed");

            // Prepare (shared) memory buffer
            // HACK: misuse error message area as array of bad file indices
            // There should be space for the max number of handled files plus the length field
            const unsigned int message_size = error_collector::singleton->in_collect_mode ? std::max((file_paths.size() / nprocs + 2) * sizeof(uint32_t), 64ul) : 0;
            mp::buf<typename mx_trait<f_type>::complex_type> buf(mx_trait<f_type>::get_complex(matrix.get()), matrix_dims, nprocs, message_size);

            auto r_time = std::chrono::high_resolution_clock::now();
            // Spawn processes if needed
            mp::run<typename mx_trait<f_type>::complex_type>(nprocs, buf, [&file_paths, &object_path, &mdims](int proc, mp::buf<typename mx_trait<f_type>::complex_type> &buf) {
                    // Read objects data into buffer
                    typename mx_trait<f_type>::complex_type *data_buf = buf.get(proc);
                    int first = buf.offset(proc);
                    int last = buf.offset(proc+1);
                    for (int obj=first; obj<last; obj++) {
                        try {
                            std::size_t offset = (obj - first) * buf.chunk_size;
                            DEBUG {
                                OUT << "Process " << proc << ": " << file_paths[obj] << '@' << offset << std::endl;
                            }
                            hdf5::file file(open_file(file_paths[obj]));
                            read_complex_matrix<f_type>(file, object_path, mdims, &data_buf[offset]);
                        } catch (std::exception &ex) {
                            error_collector::singleton->error(read_error(obj, ex.what()));
                        } catch (...) {
                            error_collector::singleton->error(read_error(obj, "mysterious read error"));
                        }
                    }
                    if (proc && error_collector::singleton->in_collect_mode) {
                        // HACK: misuse error message area as array of bad file indices for this process
                        // Layout: [ N | elem0 | elem1 | ... | elemN-1 ]
                        uint32_t *msg_area = reinterpret_cast<uint32_t*>(buf.get_msg(proc));
                        auto &bad_index = error_collector::singleton->bad_index;
                        msg_area[0] = bad_index.size();
                        std::copy(bad_index.begin(), bad_index.end(), &msg_area[1]);
                    }
                });
            // Set statistics (not correct for error collection mode)
            if (r_stat) {
                r_stat->seconds = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - r_time).count();
                r_stat->nbytes = buf.length * sizeof(typename mx_trait<f_type>::complex_type);
            }
            // Collect bad indices in error collection mode
            if (error_collector::singleton->in_collect_mode) {
                DEBUG {
                    OUT << "Collecting bad indices:" << std::endl;
                }
                for (int proc=1; proc<nprocs; proc++) {
                    // HACK: misuse error message area as array of bad file indices
                    // Layout: [ N | elem0 | elem1 | ... | elemN-1 ]
                    uint32_t *msg_area = reinterpret_cast<uint32_t*>(buf.get_msg(proc));
                    uint32_t nbad = msg_area[0];
                    DEBUG {
                        OUT << "  process " << proc << ": " << nbad << " indices" << std::endl;
                    }
                    if (nbad > (file_paths.size() / nprocs + 1))
                        throw std::runtime_error("transfer of bad indices failed");
                    for (uint32_t i=1; i<=nbad; i++)
                        error_collector::singleton->add(msg_area[i]);
                }
            }
            // Return result
            return matrix.release();
        } catch (std::exception &ex) {
            mexErrMsgIdAndTxt("ptycho:read:failed", "%s", ex.what());
            return nullptr;
        }
    }

    /*!
     * \brief Convert vector<uint32_t> to MATLAB array
     * \param vec vector to be converted
     * \return MATLAB array
     */
    mxArray* vector_to_array(std::vector<uint32_t> &vec)
    {
        mwSize vlen[2] = {1, vec.size()};
        mx_ptr<mxArray> arr(mxCreateNumericArray(2, vlen, mxUINT32_CLASS, mxREAL));
        if (! arr.get())
            throw std::runtime_error("unable to create bad indices array");
        if (! vec.empty()) {
            auto* element = mxGetUint32s(arr.get());
            if (! element)
                throw std::runtime_error("unable to retrieve pointer for bad indices array elements");
            for (auto& e : vec)
                *element++ = e+1;
        }
        return arr.release();
    }

} // namespace

namespace data_prep {

    mxArray* read_object_data(int nprocs, precision::type prec, std::array<int64_t, 2> &dims,
                              const std::string &object_path, const std::vector<std::string> &file_paths,
                              mxArray **bad_file_idx, read_stat *r_stat)
    {
        if (! file_paths.size())
            return nullptr;

        if (nprocs <= 0)
            nprocs = 1;
        if (nprocs > file_paths.size())
            nprocs = (int)file_paths.size();

        error_collector errcol(bad_file_idx != nullptr);

        mx_ptr<mxArray> result;
        {
            if (prec == precision::type::Double)
                result.reset(read_objects_parallel<double>(nprocs, dims, object_path, file_paths, r_stat));
            else
                result.reset(read_objects_parallel<float>(nprocs, dims, object_path, file_paths, r_stat));
        }

        if (errcol.in_collect_mode)
            *bad_file_idx = vector_to_array(errcol.bad_index);

        return result.release();
    }

} // namespace data_prep
