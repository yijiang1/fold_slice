#ifndef READ_DATA_THREADED_H
#define READ_DATA_THREADED_H

/*!
 * \file
 * Header file for the read Moench and Pilatus data functionality
 */

/*!
 * \brief Data preparator functions
 */
namespace data_prep {
    /*!
     * \brief read Moench or Pilatus data
     *
     * \param format data format string
     * \param nthreads number of parallel read threads
     * \param data_path paths to data files
     * \param image_size size of returned image data
     * \param roi_center region of interest center relative to measured data
     * \param prec precision of returned data
     * \return MATLAB array with resized and recentered images
     */
    mxArray* read_data_threaded(const std::string &format,
                                long nthreads,
                                const std::vector<std::string> &data_path,
                                std::vector<long> &image_size, std::vector<long> &roi_center,
                                precision::type prec);
}

#endif
