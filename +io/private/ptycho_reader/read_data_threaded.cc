/*!
 * \file
 * Support for data formats with many separate compnents that can be read using threads
 *
 * TODO: make the code safe to type inconsistencies
 */

#include <map>
#include <cerrno>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <cassert>
#include <cstdint>
#include <regex>
#include <mutex>
#include <future>
#include <atomic>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <dirent.h>
#include "tiffio.h"
#include "precision.h"
#include "mex.h"
#include "mex_helper.h"
#include "debug_helper.h"
#include "read_data_threaded.h"

namespace {
    thread_local unsigned int thread_id = 0;
    std::mutex mutex_cout;
    std::atomic<bool> thread_error;

#define LOCK(m) std::lock_guard<std::mutex> _lock(m);

    /*!
     * \brief Adapt file and destination dimensions
     *
     * \param fdims (IN) file image dimensions (2D)
     * \param mdims (IN) MATLAB destination array dimensions (> 2D) (last dimensions are image dimensions)
     * \param im_ctr (IN) image center
     * \param fstride (OUT) file image strides
     * \param mstride (OUT) destination image strides
     * \param count (OUT) elements to copy in each image dimension
     */
    void adapt_dimensions(const std::vector<long> &fdims, const std::vector<mwSize> mdims,
                          const std::vector<long> &im_ctr,
                          std::vector<long> &fstride,
                          std::vector<long> &mstride,
                          std::vector<long> &count)
    {
        const mwSize *md = &(*(mdims.end() - 2));
        const std::vector<mwSize> mh{md[0]/2, md[1]/2};
        for (unsigned int i=0; i<2; i++) {
            if (im_ctr[i] >= mh[i]) {
                mstride[i] = 0;
                fstride[i] = im_ctr[i] - mh[i];
                if (fstride[i] + md[i] > fdims[i])
                    count[i] = fdims[i] - fstride[i];
                else
                    count[i] = md[i];
            } else {
                fstride[i] = 0;
                mstride[i] = mh[i] - im_ctr[i];
                if (mstride[i] + fdims[i] <= md[i])
                    count[i] = fdims[i];
                else
                    count[i] = md[i] - mstride[i];
            }
        }
    }

    /*!
     * \brief Append element to file path
     *
     * Path separator is '/'
     *
     * \param path path appended to
     * \param elem element to append
     * \return new path with extra element
     */
    std::string path_append (const std::string &path, const std::string &elem)
    {
        std::string result(path);
        if (result.size())
            result.push_back('/');
        return result + elem;
    }

    // ----------------------------------------------------------------

    /* PILATUS CBF DATA
        _array_data.data
        ;
        --CIF-BINARY-FORMAT-SECTION--
        Content-Type: application/octet-stream;
            conversions="x-CBF_BYTE_OFFSET"
        Content-Transfer-Encoding: BINARY
        X-Binary-Size: 2499331
        X-Binary-ID: 1
        X-Binary-Element-Type: "signed 32-bit integer"
        X-Binary-Element-Byte-Order: LITTLE_ENDIAN
        Content-MD5: XoY7+gfct1+OKiJlzKOiRw==
        X-Binary-Number-of-Elements: 2476525
        X-Binary-Size-Fastest-Dimension: 1475
        X-Binary-Size-Second-Dimension: 1679
        X-Binary-Size-Padding: 4095
     */
    std::regex pilatus_ncols_regex(R"(X-Binary-Size-Fastest-Dimension: (\d+))");    //!< Regex for pilatus ncols
    std::regex pilatus_nrows_regex(R"(X-Binary-Size-Second-Dimension: (\d+))");     //!< Regex for pilatus nrows

    /*!
     * \brief Read CBF image dimensions
     *
     * \param data (IN) CBF image file data
     * \param dims (OUT) to be filled with [nrows, ncols] from the data
     */
    void cbf_dims (const std::vector<char> &data, std::vector<long> &dims)
    {
        assert(dims.size() == 2);
        {
            std::cmatch match;
            if (! std::regex_search(&data[0], &data[0] + data.size(), match, pilatus_ncols_regex))
                throw std::runtime_error("unable to find number of columns");
            dims[1] = std::stol(match[1].str());
            if (dims[1] <= 0)
                throw std::runtime_error("dimension along row is not positive");
        }
        {
            std::cmatch match;
            if (! std::regex_search(&data[0], &data[0] + data.size(), match, pilatus_nrows_regex))
                throw std::runtime_error("unable to find number of rows");
            dims[0] = std::stol(match[1].str());
            if (dims[0] <= 0)
                throw std::runtime_error("dimension along column is not positive");
        }
    }

    /*!
     * \brief Read pilatus metadata
     *
     * \param paths (IN) paths to pilatus CBF image files
     * \param dims (OUT) to be filled with [n_images, n_series(only if several series are present), nrows, ncols], the last two from the first data file
     * \param im_sz (INOUT) if empty, fill it with [nrows, ncols] from the first data file
     */
    void pilatus_read_meta(const std::vector<std::string> &paths,
                           std::vector<long> &dims,
                           std::vector<long> &im_sz)
    {
        // Take the first file to determine image dimensions size and center if not given
        if (paths.empty())
            throw std::invalid_argument("no CBF data files");
        DEBUG {
            OUT << paths.size() << " data files\n"
                << "first file " << paths[0] << std::endl;
        }
        std::vector<long> fdim{0, 0};
        {
            std::vector<char> buf(4*1024);  // 4K max header size in CBF file
            std::ifstream ifs(paths[0]);
            ifs.read(&buf[0], buf.size());
            if (ifs.fail() && !ifs.eof())
                throw std::runtime_error("unable to read first data file");
            buf.resize(ifs.gcount());
            cbf_dims(buf, fdim);
        }
        dims.resize(3);
        dims[0] = paths.size();
        dims[1] = fdim[0];
        dims[2] = fdim[1];
        if (im_sz.empty()) {
            im_sz.resize(2);
            im_sz[0] = fdim[0];
            im_sz[1] = fdim[1];
            DEBUG {
                OUT << "setting image size to " << im_sz[0] << 'x' << im_sz[1] << std::endl;
            }
        }
    }

    /*!
     * \brief Read image data from pilatus CBF file
     *
     * \param tid thread id starting from 0
     * \param destination MATLAB array image data destination
     * \param path pilatus CBF data file path
     * \param mdims MATLAB array dimensions
     * \param im_ctr image center relative to data
     * \tparam f_type float or double array element type
     */
    template<typename f_type>
    void pilatus_read_data(unsigned int tid,
                           f_type * destination,
                           const std::string &path,
                           const std::vector<mwSize> &mdims,
                           const std::vector<long> &im_ctr)
    {
        // Use Heiners method to read the data
        // Adapt dimensions for every image file
        std::vector<char> fbuf;
        {
            FILE *fin = fopen(path.c_str(), "r");
            if (! fin)
                throw std::runtime_error(std::string("unable to open file ") + path + ": " + std::strerror(errno));
            try {
                off_t fsz;
                {
                    struct stat sbuf;
                    if (fstat(fileno(fin), &sbuf) == -1) {
                        throw std::runtime_error(std::string("unable to stat file ") + path + ": " + std::strerror(errno));
                    }
                    fsz = sbuf.st_size;
                }
                fbuf.resize(fsz);
                fread(&fbuf[0], 1, fsz, fin);
                if (ferror(fin))
                    throw std::runtime_error(std::string("unable to read file ") + path + ": " + std::strerror(errno));
                fclose(fin);
            } catch (...) {
                fclose(fin);
                throw;
            }
        }
        std::vector<long> fdims{0, 0};
        cbf_dims(fbuf, fdims);
        unsigned long finger;    // compressed data index
        {
            std::vector<char> sig{ '\x0c', '\x1a', '\x04', '\xd5' };
            auto p = std::search(fbuf.begin(), fbuf.end(), sig.begin(), sig.end());
            if (p == fbuf.end())
                throw std::runtime_error(std::string("data signature not found within file " + path));
            finger = p - fbuf.begin() + 4;
        }
        unsigned long nelems = fdims[0] * fdims[1];
        std::vector<f_type> data(nelems);
        int current = 0;
        for (unsigned int i=0; i<nelems; i++) {
            if (*((uint8_t *)&fbuf[finger]) != 0x80) {    // | xx |
                current += *((int8_t *)&fbuf[finger]);
                finger += 1;
            } else if (*((uint16_t *)&fbuf[finger+1]) != 0x8000) {  // | 80 | xx | xx |
                current += *((int16_t *)&fbuf[finger+1]);
                finger += 3;
            } else {    // | 80 | 80 | 00 | xx | xx | xx | xx |
                current += *((int32_t *)&fbuf[finger+3]);
                finger += 7;
            }
            if (finger + 7 > fbuf.size())
                throw std::runtime_error(std::string("data inconsistency in file ") + path);
            //if (current < -1)   // allow value -1, which is used to mark detector gaps
            //    throw std::runtime_error(std::string("data error in file ") + path);
            data[i] = current;
        }
        fbuf.clear();
        std::vector<long> fstride{0, 0};
        std::vector<long> mstride{0, 0};
        std::vector<long> count{0, 0};
        adapt_dimensions(fdims, mdims, im_ctr, fstride, mstride, count);
        DEBUG {
            LOCK(mutex_cout);
            OUT << thread_id << ":     fstride=[" << fstride[0] << ',' << fstride[1] << "], mstride=[" << mstride[0] << ',' << mstride[1] << "], count=[" << count[0] << ',' << count[1] <<']' << std::endl;
        }
        {
            auto sz = mdims.size() - 2;
            auto msize = mdims[sz] * mdims[sz+1];
            std::memset(destination, 0, msize * sizeof(f_type));
        }
        for (unsigned long row=0; row<count[0]; row++) {
            auto col = data.begin() + fdims[1] * (row + fstride[0]) + fstride[1];
            std::copy(col, col + count[1], destination + mdims.back() * (row + mstride[0]) + mstride[1]);
        }
    }

    //-----------------------------------------------------

    /*!
     * \brief Read moench tiff image dimensions
     *
     * \param tiff_handle (IN) TIFF file handle
     * \param fdims (OUT) to be filled with [nrows, ncols]
     */
    void tiff_dims (TIFF *tiff_handle, std::vector<long> &fdims)
    {
        assert(fdims.size() == 2);
        uint32_t image_length, image_width;
        if (! TIFFGetField(tiff_handle, TIFFTAG_IMAGELENGTH, &image_length))
            throw std::runtime_error("unable to get image length");
        if (! image_length)
            throw std::runtime_error("dimension along column is not positive");
        if (! TIFFGetField(tiff_handle, TIFFTAG_IMAGEWIDTH, &image_width))
            throw std::runtime_error("unable to get image width");
        if (! image_width)
            throw std::runtime_error("dimension along row is not positive");
        fdims[0] = image_length;
        fdims[1] = image_width;
    }

    /*!
     * \brief Read moench metadata
     *
     * \param paths (IN) list of moench TIFF image data file paths
     * \param dims (OUT) data dimensions [n_images, nrows, ncols], the last two from the first data file
     * \param im_sz (INOUT) if empty, set to [nrows, ncols] from first data file
     */
    void moench_read_meta(const std::vector<std::string> &paths,
                          std::vector<long> &dims,
                          std::vector<long> &im_sz)
    {
        // Take the first file to determine image dimensions size and center if not given
        if (paths.empty())
            throw std::invalid_argument("no TIFF data files");
        DEBUG {
            OUT << paths.size() << " data files\n"
                << "first file " << paths[0] << std::endl;
        }
        std::vector<long> fdim(2);
        {
            TIFF *tiff_handle = TIFFOpen(paths[0].c_str(), "r");
            if (! tiff_handle)
                throw std::runtime_error("unable to open first image data file");
            try {
                tiff_dims(tiff_handle, fdim);
            } catch (...) {
                TIFFClose(tiff_handle);
                throw;
            }
            TIFFClose(tiff_handle);
        }
        dims.resize(3);
        dims[0] = paths.size();
        dims[1] = fdim[0];
        dims[2] = fdim[1];
        if (im_sz.empty()) {
            im_sz.resize(2);
            im_sz[0] = fdim[0];
            im_sz[1] = fdim[1];
            DEBUG {
                OUT << "setting image size to " << im_sz[0] << 'x' << im_sz[1] << std::endl;
            }
        }
    }

    /*!
     * \brief Read TIFF image data
     *
     * \param tf TIFF file descriptor
     * \param buf char buffer with enough space for a data strip
     * \param image image data buffer
     * \param num_strip number of data strips
     * \tparam sample_type image data sample type
     * \tparam result_type image buffer data type
     */    
    template<typename sample_type, typename result_type>
    void tiff_read(TIFF *tf, std::vector<char> &buf, std::vector<result_type> &image, tstrip_t num_strips)
    {
        std::uint32_t idx = 0;
        for (tstrip_t strip=0; strip<num_strips; strip++) {
            tsize_t nbytes = TIFFReadEncodedStrip(tf, strip, buf.data(), buf.size());
            if (nbytes < 0)
                throw std::runtime_error("unable to read strip from tiff file");
            for (tsize_t i=0; i<nbytes; i+=sizeof(sample_type))
                image[idx++] = *((sample_type *)&buf[i]);
        }
        if (idx != image.size())
            throw std::runtime_error("tiff image data size mismatch");
    }


    /*!
     * \brief Read moench TIFF image data into MATLAB array
     *
     * \param tid thread id, starting from 0
     * \param destination MATLAB array destination data buffer
     * \param path tiff file path
     * \param mdims MTLAB array dimensions
     * \param im_ctr image center relative to data
     */
    template<typename f_type>
    void moench_read_data (unsigned int tid,
                           f_type * destination,
                           const std::string &path,
                           const std::vector<mwSize> &mdims,
                           const std::vector<long> &im_ctr)
    {
        /*
        std::vector<long> fdims{0, 0};
        std::vector<float> fbuf(0);
        {
            TIFF *tiff_handle = TIFFOpen(path.c_str(), "r");
            if (! tiff_handle)
                throw std::runtime_error(std::string("unable to open file ") + path);
            try {
                tiff_dims(tiff_handle, fdims);
                long nelems = fdims[0] * fdims[1];
                fbuf.resize(nelems);
                tmsize_t res, sz = fdims[1] * sizeof(float);
                for (uint32_t strip=0; strip<fdims[0]; strip++) {
                    float *pos = &fbuf[fdims[1] * strip];
                    res = TIFFReadRawStrip(tiff_handle, strip, pos, sz);
                    if (res != sz)
                        throw std::runtime_error("unable to read tiff data");
                }
            } catch (...) {
                TIFFClose(tiff_handle);
                throw;
            }
            TIFFClose(tiff_handle);
        }
        */
        std::vector<long> fdims{0, 0};
        std::vector<f_type> fbuf(0);
        {
            TIFF *tiff_handle = TIFFOpen(path.c_str(), "r");
            if (! tiff_handle)
                throw std::runtime_error(std::string("unable to open file ") + path);
            try {
                tiff_dims(tiff_handle, fdims);
                long nelems = fdims[0] * fdims[1];
                fbuf.resize(nelems);
                std::uint16_t bps;
                if (TIFFGetField(tiff_handle, TIFFTAG_BITSPERSAMPLE, &bps) != 1)
                    throw std::runtime_error(std::string("unable to read number of bits per sample for file ") + path);
                std::uint16_t format;
                if (TIFFGetFieldDefaulted(tiff_handle, TIFFTAG_SAMPLEFORMAT, &format) != 1)
                    throw std::runtime_error(std::string("unable to read sample format for file") + path);
                tstrip_t num_strips = TIFFNumberOfStrips(tiff_handle);
                tsize_t strip_sz = TIFFStripSize(tiff_handle);
                std::vector<char> cbuf(strip_sz);
                switch (format) {
                case 1:
                    switch (bps) {
                    case 8:
                        tiff_read<std::uint8_t, f_type>(tiff_handle, cbuf, fbuf, num_strips); break;
                    case 16:
                        tiff_read<std::uint16_t, f_type>(tiff_handle, cbuf, fbuf, num_strips); break;
                    case 32:
                        tiff_read<std::uint32_t, f_type>(tiff_handle, cbuf, fbuf, num_strips); break;
                    default:
                        throw std::runtime_error(std::string("unsupported number of bits per unsigned integer sample in file ") + path);
                    }
                    break;
                case 2:
                    switch (bps) {
                    case 8:
                        tiff_read<std::int8_t, f_type>(tiff_handle, cbuf, fbuf, num_strips); break;
                    case 16:
                        tiff_read<std::int16_t, f_type>(tiff_handle, cbuf, fbuf, num_strips); break;
                    case 32:
                        tiff_read<std::int32_t, f_type>(tiff_handle, cbuf, fbuf, num_strips); break;
                    default:
                        throw std::runtime_error(std::string("unsupported number of bits per integer sample in file ") + path);
                    }
                    break;
                case 3:
                    switch (bps) {
                    case 32:
                        tiff_read<float, f_type>(tiff_handle, cbuf, fbuf, num_strips); break;
                    case 64:
                        tiff_read<double, f_type>(tiff_handle, cbuf, fbuf, num_strips); break;
                    default:
                        throw std::runtime_error(std::string("unsupported number of bits per ieee sample in file ") + path);
                    }
                    break;
                default:
                    throw std::runtime_error(std::string("unsupported sample format in file ") + path);
                }
            } catch (...) {
                TIFFClose(tiff_handle);
                throw;
            }
            TIFFClose(tiff_handle);
        }
        std::vector<long> fstride{0, 0};
        std::vector<long> mstride{0, 0};
        std::vector<long> count{0, 0};
        adapt_dimensions(fdims, mdims, im_ctr, fstride, mstride, count);
        DEBUG {
            LOCK(mutex_cout);
            OUT << thread_id << ":     fstride=[" << fstride[0] << ',' << fstride[1] << "], mstride=[" << mstride[0] << ',' << mstride[1] << "], count=[" << count[0] << ',' << count[1] <<']' << std::endl;
        }
        {
            unsigned long msize = *(mdims.end() - 2) * mdims.back();
            std::memset(destination, 0, msize * sizeof(f_type));
        }
        for (unsigned long row=0; row<count[0]; row++) {
            auto col = fbuf.begin() + fdims[1] * (row + fstride[0]) + fstride[1];
            std::copy(col, col + count[1], destination + mdims.back() * (row + mstride[0]) + mstride[1]);
        }
    }

    //-----------------------------------------------------

    /*!
     * \brief Read data files in parallel
     *
     * \param format data format name
     * \param nthreads desired number of threads
     * \param paths paths to detector data files
     * \param image_size desired image size [nrows, ncols]
     * \param roi_center desired image center [row, col] relative to data
     * \return MATLAB array [n_images, n_series(only if several are present), nrows, ncols]
     */
    template<typename f_type>
    mxArray* read_data_parallel(const std::string &format,
                                long nthreads,
                                const std::vector<std::string> &paths,
                                std::vector<long> &image_size, std::vector<long> &roi_center)
    {
        struct read_func final {
            void (*read_meta) (const std::vector<std::string> &paths,
                               std::vector<long> &dims,
                               std::vector<long> &im_sz);
            void (*read_data) (unsigned int tid,
                               f_type * destination,
                               const std::string &path,
                               const std::vector<mwSize> &mdims,
                               const std::vector<long> &im_ctr);
        } detector_functions[2] = {
            { pilatus_read_meta, pilatus_read_data<f_type> },   // pilatus
            { moench_read_meta, moench_read_data<f_type> }      // moench
        };

        const std::map<std::string, read_func&> dfunc {
            { "pilatus",  detector_functions[0] },
            { "moench", detector_functions[1] }
        };

        const auto elem = dfunc.find(format);
        if (elem == dfunc.end())
            throw std::invalid_argument("unsupported data format");
        const auto &func = elem->second;
        std::vector<long> fdims;
        func.read_meta(paths, fdims, image_size);
        if (fdims.size() < 3)
            throw std::invalid_argument("inconsistent data (dimensionality too low)");
        if (paths.size() != fdims[0])
            throw std::invalid_argument("inconsistent data (number of paths / array dimension mismatch)");
        if (image_size.size() != 2)
            throw std::invalid_argument("bad image size, must be two dimensional");
        if (roi_center.empty()) {
            roi_center.resize(2);
            roi_center[0] = fdims[fdims.size() - 2] / 2;
            roi_center[1] = fdims.back() / 2;
        }
        std::vector<mwSize> mdims(fdims.size());
        std::transform(fdims.begin(), fdims.end(), mdims.begin(), [](long e)->mwSize { return (mwSize)e; });
        auto sz = mdims.size() - 2;
        mdims[sz] = image_size[0];
        mdims[sz + 1] = image_size[1];
        DEBUG {
            OUT << "fdims  " << fdims[0] << 'x' << fdims[1] << 'x' << fdims[2] << '\n'
                << "mdims  " << mdims[0] << 'x' << mdims[1] << 'x' << mdims[2] << '\n'
                << "center " << roi_center[0] << 'x' << roi_center[1] << std::endl;
        }
        mx_ptr<mxArray> matrix;
        {
            std::vector<mwSize> rdims(mdims);
            std::reverse(rdims.begin(), rdims.end());
            matrix.reset(mxCreateNumericArray(rdims.size(), &rdims[0], mx_trait<f_type>::class_id, mxREAL));
        }
        if (! matrix.get())
            throw std::runtime_error("matrix creation failed");
        f_type *destination = mx_trait<f_type>::get(matrix.get());

        if (nthreads > paths.size())
            nthreads = paths.size();
        std::vector<std::future<void>> threads;
        for (unsigned int tid=0; tid<nthreads; tid++) {
            threads.push_back(std::async(std::launch::async, [nthreads, tid, &destination, &func, &paths, &fdims, &mdims, &roi_center]() {
                thread_id = tid;
                const unsigned int npaths = (paths.size() + nthreads - 1) / nthreads;
                const unsigned int first_full = (nthreads * npaths) - paths.size();
                const unsigned int first = (tid * npaths) - (tid < first_full ? tid : first_full);
                const unsigned int last = first + npaths - (tid < first_full ? 1 : 0);
                const auto pdim = &mdims[mdims.size() - 2];
                const auto sz = pdim[0] * pdim[1];
                DEBUG {
                    LOCK(mutex_cout);
                    OUT << thread_id << ": range " << first << '-' << last << ", sz=" << sz << std::endl;
                }
                for (unsigned int i=first; i<last; i++) {
                    if (thread_error.load())
                        throw std::runtime_error("received stop signal");
                    DEBUG {
                        LOCK(mutex_cout);
                        OUT << thread_id << ": " << paths[i] << std::endl;
                    }
                    func.read_data(tid, &destination[i * sz], paths[i], mdims, roi_center);
                }
            }));
        }

        bool success{true};
        for (unsigned int tid=0; tid<nthreads; tid++) {
            try {
                threads[tid].wait();
                threads[tid].get();
            } catch (const std::exception &ex) {
                mexPrintf("Thread %u error: %s\n", tid, ex.what());
                thread_error.store(false);
            }
        }

        if (thread_error)
            throw std::runtime_error("Failed to read data!");

        return matrix.release();
    }

}

namespace data_prep {
    mxArray* read_data_threaded(const std::string &format,
                                long nthreads,
                                const std::vector<std::string> &data_path,
                                std::vector<long> &image_size, std::vector<long> &roi_center,
                                precision::type prec)
    {
        if (prec == precision::type::Double)
            return read_data_parallel<double>(format, nthreads, data_path, image_size, roi_center);
        else
            return read_data_parallel<float>(format, nthreads, data_path, image_size, roi_center);
    }
}
