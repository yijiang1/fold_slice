#ifndef READ_OBJECT_DATA_H
#define READ_OBJECT_DATA_H

/*!
 * \file
 * Declaration of the read_object_data function
 */

/*!
 * \brief Data preparator functions
 */
namespace data_prep {

    /*!
     * \brief Statistics collection
     */
    struct read_stat final {
        double seconds;         //!< read access time in seconds
        unsigned long nbytes;   //!< number of bytes transferred
    };

    /*!
     * \brief Read object data
     *
     * The function will either create the bad indices vector or return with exception and error message.
     *
     * \param nprocs Number of read processes
     * \param prec Single or Double precision
     * \param dims Object size (nrows, ncols)
     * \param object_path Path to object within HDF5 file
     * \param file_paths Paths to HDF5 files, each containing a diffraction pattern at object_path
     * \param bad_file_idx Indices of nonreadable files array, to be created if not NULL
     * \param r_stat Fill with info if not null
     * \return MATLAB array with object data
     */
    mxArray* read_object_data(int nprocs, precision::type prec, std::array<int64_t, 2> &dims,
                              const std::string &object_path, const std::vector<std::string> &file_paths,
                              mxArray **bad_file_idx,
                              read_stat *r_stat = nullptr);

} // namespace data_prep

#endif
