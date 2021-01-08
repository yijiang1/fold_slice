#ifndef READ_EIGER_DATA_H
#define READ_EIGER_DATA_H

/*!
 * \file
 * Header file for the read Eiger data functionality
 */

/*!
 * \brief Data preparator functions
 */
namespace data_prep {

    /*!
     * \brief read Eiger data
     *
     * \param nprocs number of parallel read processes
     * \param data_path paths to Eiger data files
     * \param data_location location within file, either 1 for all paths or one for each path
     * \param image_size size of returned image data
     * \param roi_center region of interest center relative to measured data
     * \param prec precision of returned data
     * \return MATLAB array with resized and recentered images
     * \pre The 1:N or N:N relationship between data_location and data_path must hold
     */
    mxArray* read_eiger_data(long nprocs,
                             const std::vector<std::string> &data_path,
                             const std::vector<std::string> &data_location,
                             std::vector<long> &image_size, std::vector<long> &roi_center,
                             precision::type prec);
}

#endif
