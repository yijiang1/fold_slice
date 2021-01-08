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
#include <cstdio>
#include <functional>
#include <algorithm>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/wait.h>
#include "hdf5_helper.h"
#include "mex.h"
#include "debug_helper.h"
#include "precision.h"
#include "multi_processing.h"

namespace {

    #include "mex_helper.h"

    /*!
     * \defgroup dataLayoutKind Kind of data layout
     * @{
     */
    constexpr int kind_11 = 0;                   //!< 1 file - 1 dataset   (new eiger)
    constexpr int kind_N1 = 1;                   //!< N files - 1 dataset  (old eiger)
    constexpr int kind_1N = 2;                   //!< 1 file - N datasets  (ESRF)
    constexpr int kind_NN = 3;                   //!< N (file, dataset) pairs
    /*! @} */

    /*!
     * \brief Vector to string
     * Outputs a string containing v0[,v1[,v2...]]
     * \param v vector
     * \param separator separator between vector elements, defaults to comma
     * \tparam T vector element type
     */
    template <typename T>
    std::string vec_to_str(const std::vector<T> &v, char separator=',')
    {
        std::ostringstream oss;
        oss << v[0];
        for (std::size_t i=1; i<v.size(); i++)
            oss << separator << v[i];
        return oss.str();
    }

    /*!
     * \brief Adapt data space layout
     *
     * Adapt the memory and file dataset layouts by selecting apropriate hyperslabs
     *
     * \param mspace memory data space layout (must fit to parameters first, last)
     * \param fspace file dataset layout
     * \param mdim memory space dimensions (must fit to parameters first last if kind=kind_11)
     * \param fdim file dataset dimensions
     * \param first first chunk of file space to include
     * \param last one past last chunk of file space to include
     * \param roi_center image center
     * \tparam kind kind of data layout
     */
    template<int kind>
    void adapt_spaces(hdf5::dataspace &mspace, hdf5::dataspace &fspace,
                      const std::vector<hsize_t> &mdim, const std::vector<hsize_t> &fdim,
                      int first , int last,
                      const std::vector<long> &roi_center)
    {
        if (mdim.size() <= roi_center.size())
            throw std::invalid_argument("result matrix must have more dimensions than the image");

        std::vector<hsize_t> fstart(fdim.size());
        std::vector<hsize_t> mstart(fdim.size());
        std::vector<hsize_t> fcount(fdim);
        std::vector<hsize_t> mcount(mdim);
        bool mselect = false;
        bool fselect = false;

        if (kind == kind_11) {
            // adapt first dimension (chunk dimension)
            if ((first != 0) || (last != fdim[0])) {
                fstart[0] = first;
                fcount[0] = last - first;
                fselect = true;
                // DEBUG {
                //     OUT << "  adapted chunk dimension: " << fstart[0] << '+' << fcount[0] << std::endl;
                // }
            }
        }

        // adapt last dimensions (normally detector image dimensions)
        {
            auto fi = fdim.size();
            auto mi = mdim.size();
            auto i = roi_center.size();
            do {
                i--; fi--; mi--;
                hsize_t fh = (fdim[fi] + 1) / 2;
                hsize_t mh = (mdim[mi] + 1) / 2;
                hsize_t rc = roi_center[i];
                if (rc<0 || rc>=fdim[fi])
                    throw std::invalid_argument("roi_center out of image bounds");
                if (mh <= rc) {
                    fstart[fi] = rc - mh;
                    if (fstart[fi] + mdim[mi] >= fdim[fi]) {
                        fcount[fi] = fdim[fi] - fstart[fi];
                        mcount[mi] = fcount[fi];
                        mselect = mselect || (mcount[mi] != mdim[mi]);
                    } else {
                        fcount[fi] = mdim[mi];
                    }
                    fselect = fselect || ((fcount[fi] != fdim[fi]) || fstart[fi]);
                } else {
                    if (rc + mh <= fdim[fi]) {
                        mcount[mi] = fcount[fi] = rc + mh;
                        mstart[mi] = mdim[mi] - mcount[mi];
                        fselect = fselect || (fcount[fi] != fdim[fi]);
                    } else {
                        mcount[mi] = fcount[fi];
                        mstart[mi] = rc - mh;
                    }
                    mselect = true;
                }
                // DEBUG {
                //     OUT << "  adapting image dimension " << i << ": fdim=" << fdim[fi] << '/' << fh << ", mdim=" << mdim[mi] << '/' << mh << ", rc=" << rc << "\n  fsel=" << fselect << ", msel=" << mselect << std::endl;
                // }
            } while(i);
        }
        if (mselect) {
            DEBUG {
                OUT << "memory subspace start=" << vec_to_str(mstart) << " count=" << vec_to_str(mcount) << "\nfile space is " << vec_to_str(fdim) << std::endl;
            }
            if (H5Sselect_hyperslab(mspace, H5S_SELECT_SET, &mstart[0], nullptr, &mcount[0], nullptr) < 0)
                throw hdf5::exception("unable to set hdf5 memory space");
        }
        if (fselect) {
            DEBUG {
                OUT << "file subspace start=" << vec_to_str(fstart) << " count=" << vec_to_str(fcount) << "\nmem space is " << vec_to_str(mdim) << std::endl;
            }
            if (H5Sselect_hyperslab(fspace, H5S_SELECT_SET, &fstart[0], nullptr, &fcount[0], nullptr) < 0)
                throw hdf5::exception("unable to set hdf5 dataset space");
        }
    }

    /*!
     * \brief Read object from file into buffer
     * \param file HDF5 file
     * \param dataset_path path to dataset
     * \param first first buffer chunk to read
     * \param last last buffer chunk to read
     * \param mdims buffer dimensions
     * \param ddims dataset dimensions
     * \param buf pointer to parent acessible buffer
     * \param roi_center image center
     * \tparam f_type buffer element type
     * \tparam kind kind of data layout
     */
    template <typename f_type, int kind>
    void read_data_chunk(hdf5::file &file, const std::string &dataset_path,
                         int first, int last,
                         const std::vector<mwSize> &mdims,
                         const std::vector<hsize_t> &ddims,
                         f_type *buf,
                         const std::vector<long> &roi_center)
    {
        // Open dataset and data space
        hdf5::dataset dataset(H5Dopen(file, dataset_path.c_str(), H5P_DEFAULT));
        if (! dataset.valid())
            throw hdf5::exception("unable to open dataset");
        hdf5::dataspace dataspace(H5Dget_space(dataset));
        if (! dataspace.valid())
            throw hdf5::exception("unable to open data space for dataset");

        // Read and check data space
        std::vector<hsize_t> sdims(ddims.size()); // source dimensions
        {
            int ndims = H5Sget_simple_extent_ndims(dataspace);
            if (ndims != ddims.size())
                throw hdf5::exception("data space dimension mismatch");
            if (H5Sget_simple_extent_dims(dataspace, &sdims[0], nullptr) != ndims)
                throw hdf5::exception("unable to get dimension sizes");
            for (unsigned int i=0; i<ddims.size()-2; i++)
                if (sdims[i] != ddims[i])
                    throw std::runtime_error("data space dimension inconsistency");
        }

        std::vector<hsize_t> tdims(sdims);  // target dimensions
        if (kind == kind_11) {
            auto i = sdims.size() - mdims.size();
            decltype(i) j = 0;
            for (; i<sdims.size(); i++, j++)
                tdims[i] = mdims[j];
            tdims[0] = last - first;
        } else {
            auto i = sdims.size();
            auto j = mdims.size();
            while ((i > 0) && (j > 1)) { // don't copy first matrix dimension (number of files/datasets)
                i--; j--;
                tdims[i] = mdims[j];
            }
            while (i) { // check that remaining data file dimensions are 1
                i--;
                if (tdims[i] != 1)
                    throw std::runtime_error("data space dimension error");
            }
        }
        DEBUG {
            OUT << "init mem space: " << vec_to_str(tdims, 'x') << " at address " << buf << std::endl;
        }
        hdf5::dataspace memspace(H5Screate_simple(tdims.size(), &tdims[0], nullptr));
        if (! memspace.valid())
            throw hdf5::exception("unable to create memory data space");
        adapt_spaces<kind>(memspace, dataspace, tdims, sdims, first, last, roi_center);
        // Read dataset
        if (H5Dread(dataset, hdf5::type_trait<f_type>::type, memspace, dataspace, H5P_DEFAULT, buf) < 0) {
            DEBUG {
                std::FILE *out = std::fopen(std::getenv("PTYCHO_READ_DEBUG"), "a+");
                if (! out) {
                    OUT << "unable to open file stream: " << std::strerror(errno) << std::endl;
                } else {
                    H5Eprint(H5E_DEFAULT, out);
                    std::fclose(out);
                }
            }
            throw hdf5::exception("unable to read dataset");
        }

        DEBUG {
            OUT << "read finished - [0]=" << buf[0] << std::endl;
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
            throw hdf5::exception(std::string("unable to open file ") + file_path);
        return file.grab();
    }

    /*!
     * \brief Drop leading dimensions of size 1
     * \param dims vector for which leading dimensions of size 1 beyond the second dimension will be dropped
     * \tparam T type of dimension size
     */
    template <typename T>
    void collapse_dims (std::vector<T> &dims) noexcept
    {
        int i=0;
        for (; i<dims.size()-2; i++) {
            if (dims[i] != 1)
                break;
        }
        int j=0;
        for (; i<dims.size(); i++, j++)
            dims[j] = dims[i];
        return dims.resize(j);
    }

    /*!
     * \brief Get HDF5 data size
     *
     * \param file HDF5 file
     * \param dataset_path dataset path
     * \param dims filled with dataset dimensions
     */
    void get_data_size(hdf5::file &file, const std::string &dataset_path, std::vector<hsize_t> &dims)
    {
        hdf5::dataset dataset(H5Dopen(file, dataset_path.c_str(), H5P_DEFAULT));
        if (! dataset.valid())
            throw hdf5::exception("unable to open dataset");
        hdf5::dataspace dataspace(H5Dget_space(dataset));
        if (! dataspace.valid())
            throw hdf5::exception("unable to get data space");
        int ndims = H5Sget_simple_extent_ndims(dataspace);
        if (ndims <= 0)
            throw hdf5::exception("wrong number of dimensions in dataset");
        dims.resize(ndims);
        if (H5Sget_simple_extent_dims(dataspace, &dims[0], nullptr) != ndims)
            throw hdf5::exception("unable to get dimension sizes");
    }

    /*!
     * \brief Traits for kind of data
     * \tparam kind kind of data layout
     */
    template<int kind>
    struct ktrait final {};

    /*!
     * \brief Traits for new eiger data
     */
    template<>
    struct ktrait<kind_11> final {
        /*!
         * \brief Assign matrix dimensions from data dimensions
         * \param mdims (OUT) matrix dimensions to be assigned
         * \param ddims (IN) data dimensions
         * \param data_path (IN) path of data file
         * \param data_location (IN) location of dataset
         */
        static void assign_matrix_dims(std::vector<mwSize> &mdims, const std::vector<hsize_t> &ddims, const std::vector<std::string> &data_path, const std::vector<std::string> &data_location)
        {
            auto ndims = ddims.size();
            mdims.resize(ndims);
            for (decltype(ndims) i=0; i<ndims; i++) {
                mdims[i] = static_cast<mwSize>(ddims[i]);
            }
        }

        /*!
         * \brief Path range start offset
         * \param first start offset in work range
         * \param last end offset in work range (exclusive)
         * \return start path start range offset
         */
        static int start_path(int first, int last)
        {
            return 0;
        }

        /*!
         * \brief Path range end offset
         * \param first start offset in work range
         * \param last end offset in work range (exclusive)
         * \return path range end offset (exclusive)
         */
        static int end_path(int first, int last)
        {
            return 1;
        }

        /*!
         * \brief Location range start offset
         * \param f file index
         * \param first start offset in work range
         * \param last end offset in work range (exclusive)
         * \return start location start range offset
         */
        static int start_location(int f, int first, int last)
        {
            return 0;
        }

        /*!
         * \brief Location range end offset
         * \param f file index
         * \param first start offset in work range
         * \param last end offset in work range (exclusive)
         * \return start location end range offset
         */
        static int end_location(int f, int first, int last)
        {
            return 1;
        }

        /*!
         * \brief Buffer chunk offset
         * \param f file index
         * \param d dataset index
         * \param first start offset in work range
         * \return offset into chunk buffer
         */
        static std::size_t buf_offset(int f, int d, int first)
        {
            return 0;
        }
    };

    /*!
     * \brief Traits for old eiger data
     */
    template<>
    struct ktrait<kind_N1> final {
        static void assign_matrix_dims(std::vector<mwSize> &mdims, const std::vector<hsize_t> &ddims, const std::vector<std::string> &data_path, const std::vector<std::string> &data_location)
        {
            auto ndims = ddims.size();
            mdims.resize(ndims+1);
            mdims[0] = data_path.size();
            for (decltype(ndims) i=0; i<ndims; i++) {
                mdims[i+1] = static_cast<mwSize>(ddims[i]);
            }
        }

        static int start_path(int first, int last)
        {
            return first;
        }

        static int end_path(int first, int last)
        {
            return last;
        }

        static int start_location(int f, int first, int last)
        {
            return 0;
        }

        static int end_location(int f, int first, int last)
        {
            return 1;
        }

        static std::size_t buf_offset(int f, int d, int first)
        {
            return f - first;
        }
    };

    /*!
     * \brief Traits for ESRF data
     */
    template<>
    struct ktrait<kind_1N> final {
        static void assign_matrix_dims(std::vector<mwSize> &mdims, const std::vector<hsize_t> &ddims, const std::vector<std::string> &data_path, const std::vector<std::string> &data_location)
        {
            auto ndims = ddims.size();
            mdims.resize(ndims+1);
            mdims[0] = data_location.size();
            for (decltype(ndims) i=0; i<ndims; i++) {
                mdims[i+1] = static_cast<mwSize>(ddims[i]);
            }
        }

        static int start_path(int first, int last)
        {
            return 0;
        }

        static int end_path(int first, int last)
        {
            return 1;
        }

        static int start_location(int f, int first, int last)
        {
            return first;
        }

        static int end_location(int f, int first, int last)
        {
            return last;
        }

        static std::size_t buf_offset(int f, int d, int first)
        {
            return d - first;
        }
    };

    /*!
     * \brief Traits for N (file, dataset) pairs
     */
    template<>
    struct ktrait<kind_NN> final {
        static void assign_matrix_dims(std::vector<mwSize> &mdims, const std::vector<hsize_t> &ddims, const std::vector<std::string> &data_path, const std::vector<std::string> &data_location)
        {
            auto ndims = ddims.size();
            mdims.resize(ndims+1);
            mdims[0] = data_path.size();
            for (decltype(ndims) i=0; i<ndims; i++) {
                mdims[i+1] = static_cast<mwSize>(ddims[i]);
            }
        }

        static int start_path(int first, int last)
        {
            return first;
        }

        static int end_path(int first, int last)
        {
            return last;
        }

        static int start_location(int f, int first, int last)
        {
            return f;
        }

        static int end_location(int f, int first, int last)
        {
            return f+1;
        }

        static std::size_t buf_offset(int f, int d, int first)
        {
            return f - first;
        }
    };

    /*!
     * \brief Read in objects with multiple processes
     *
     * \param nprocs number of read processes
     * \param data_path paths to hdf5 files
     * \param data_location paths to datasets within hdf5 files
     * \param image_size two dimensional image size
     * \param roi_center center of two dimensional image
     * \tparam f_type float or double result type
     * \tparam kind kind of data layout (see above)
     * \return MATLAB result array
     */
    template <typename f_type, int kind>
    mxArray* read_data_parallel(int nprocs,
                                const std::vector<std::string> &data_path,
                                const std::vector<std::string> &data_location,
                                std::vector<long> &image_size, std::vector<long> &roi_center)
    {
        using kt = ktrait<kind>;
        try {
            std::vector<mwSize> mdims;      // matrix dimensions
            std::vector<hsize_t> ddims;     // dataset dimensions

            // Assign matrix dimension values
            {
                hdf5::file file(open_file(data_path[0]));
                if (! file.valid())
                    throw hdf5::exception("unable to open eiger data file");
                get_data_size(file, data_location[0], ddims);
                auto ndims = ddims.size();
                if (ndims < 2)
                    throw std::invalid_argument("dataset must have at least two dimensions");
                DEBUG {
                    OUT << "data size: " << vec_to_str(ddims, 'x') << std::endl;
                }

                auto sz = ndims - 2;
                if (image_size.empty()) {
                    image_size.resize(2);
                    image_size[0] = ddims[sz];
                    image_size[1] = ddims[sz+1];
                    DEBUG {
                        OUT << "setting image size to " << image_size[0] << 'x' << image_size[1] << std::endl;
                    }
                } else if (image_size.size() != 2)
                    throw std::invalid_argument("image_size must be two dimensional");

                if (roi_center.empty()) {
                    roi_center.resize(2);
                    roi_center[0] = ddims[sz] / 2;
                    roi_center[1] = ddims[sz+1] / 2;
                    DEBUG {
                        OUT << "setting roi center to " << roi_center[0] << 'x' << roi_center[1] << std::endl;
                    }
                } else if (roi_center.size() != 2)
                    throw std::invalid_argument("roi_center must be two dimensional");
            }

            kt::assign_matrix_dims(mdims, ddims, data_path, data_location);
            {
                auto ndims = mdims.size() - 2;
                mdims[ndims] = image_size[0];
                mdims[ndims+1] = image_size[1];
            }
            collapse_dims(mdims);
            {
                mwSize sz = 1;
                for (const auto &d : mdims)
                    sz *= d;
                if (sz <= 0)
                    throw std::invalid_argument("empty result matrix");
            }

            // Create space for MATLAB array
            DEBUG {
                OUT << "creating MATLAB array: " << vec_to_str(mdims, 'x') << std::endl;
            }
            mx_ptr<mxArray> matrix;
            {
                std::vector<mwSize> rdims(mdims);
                std::reverse(rdims.begin(), rdims.end());
                matrix.reset(mxCreateNumericArray(rdims.size(), &rdims[0], mx_trait<f_type>::class_id, mxREAL));
            }
            if (! matrix.get())
                throw std::runtime_error("matrix creation failed");

            // Prepare (shared) memory buffer
            if (mdims[0] < nprocs)
                nprocs = mdims[0];
            mp::buf<f_type> buf(mx_trait<f_type>::get(matrix.get()), mdims, nprocs);

//             auto r_time = std::chrono::high_resolution_clock::now();
            // Read objects data into buffer
            mp::run<f_type>(nprocs, buf, [&data_path, &data_location, &mdims, &ddims, &roi_center](int proc, mp::buf<f_type> &buf) {
                    int first = buf.offset(proc);
                    int last = buf.offset(proc+1);
                    DEBUG {
                        OUT << "Process " << proc << ": " << first << ".." << last << std::endl;
                    }
                    for (int f=kt::start_path(first, last); f<kt::end_path(first, last); f++) {
                        for (int d=kt::start_location(f, first, last); d<kt::end_location(f, first, last); d++) {
                            hdf5::file file(open_file(data_path[f]));
                            f_type *data_buf = buf.get(proc) + kt::buf_offset(f, d, first) * buf.chunk_size;
                            read_data_chunk<f_type, kind>(file, data_location[d], first, last, mdims, ddims, data_buf, roi_center);
                        } // dataset locations
                    } // file paths
                });
//             // Set statistics
//             if (r_stat) {
//                 r_stat->seconds = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - r_time).count();
//                 r_stat->nbytes = buf.length * sizeof(typename mx_trait<f_type>::complex_type);
//             }
            // Return result
            return matrix.release();
        } catch (std::exception &ex) {
            mexErrMsgIdAndTxt("ptycho:read:failed", "%s", ex.what());
            return nullptr;
        }
    }

    /*!
     * \brief Read in objects with multiple processes
     *
     * \param nprocs number of read processes
     * \param data_path paths to hdf5 files
     * \param data_location paths to datasets within hdf5 files
     * \param image_size two dimensional image size
     * \param roi_center center of two dimensional image
     * \tparam f_type float or double result type
     * \return MATLAB result array
     */
    template <typename f_type>
    mxArray* read_data_parallel(int nprocs,
                                const std::vector<std::string> &data_path,
                                const std::vector<std::string> &data_location,
                                std::vector<long> &image_size, std::vector<long> &roi_center)
    {
        if (data_path.size() == 1) {
            if (data_location.size() == 1) {
                DEBUG {
                    OUT << "11: new eiger" << std::endl;
                }
                return read_data_parallel<f_type, kind_11>(nprocs, data_path, data_location, image_size, roi_center);
            } else {
                DEBUG {
                    OUT << "1N: esrf" << std::endl;
                }
                return read_data_parallel<f_type, kind_1N>(nprocs, data_path, data_location, image_size, roi_center);
            }
        } else {
            if (data_location.size() == 1) {
                DEBUG {
                    OUT << "N1: old eiger" << std::endl;
                }
                return read_data_parallel<f_type, kind_N1>(nprocs, data_path, data_location, image_size, roi_center);
            } else {
                DEBUG {
                    OUT << "NN: file/dataset pairs" << std::endl;
                }
                return read_data_parallel<f_type, kind_NN>(nprocs, data_path, data_location, image_size, roi_center);
            }
        }
    }

} // namespace

namespace data_prep {

    mxArray* read_eiger_data(long nprocs,
                             const std::vector<std::string> &data_path,
                             const std::vector<std::string> &data_location,
                             std::vector<long> &image_size, std::vector<long> &roi_center,
                             precision::type prec)
    {
        if (nprocs <= 0)
            nprocs = 1;

        mx_ptr<mxArray> result;
        {
            if (prec == precision::type::Double)
                result.reset(read_data_parallel<double>(nprocs, data_path, data_location, image_size, roi_center));
            else
                result.reset(read_data_parallel<float>(nprocs, data_path, data_location, image_size, roi_center));
        }
        return result.release();
    }

} // namespace data_prep
