/*!
 * \file
 * Read Measurement Data
 *
 * This file contains the main mex code for reading measurement data into MATLAB memory.
 */

#include <string>
#include <memory>
#include <vector>
#include <stdexcept>
#include <fstream>
#include "mex.h"
#include "precision.h"
#include "read_eiger_data.h"
#include "read_data_threaded.h"
#include "debug_helper.h"

namespace {

    namespace p {
        constexpr char format[] = "extension";                          //!< extension field
        constexpr char image_size[] = "asize";                          //!< image size field
        constexpr char roi_center[] = "ctr";                            //!< ROI center field
        constexpr char scan_number[] = "scan_number";                   //!< scan number field
        constexpr char data_path[] = "data_path";                       //!< data path field
        constexpr char data_location[] = "data_location";               //!< data location field
        constexpr char data_prefix[] = "data_prefix";                   //!< data prefix field
        constexpr char scan_string_format[] = "scan_string_format";     //!< scan string format field
        constexpr char precision[] = "precision";                       //!< precision string ('single' or 'double')
        constexpr char nthreads[] = "nthreads";                         //!< number of parallel threads or processes
    } // namespace p

    #include "env_helper.h"
    #include "mex_helper.h"

    /*!
     * \brief Get field from struct
     *
     * \param pa Pointer to MATLAB struct array
     * \param field_name Field name
     * \return Pointer to MATALAB array of field value
     */
    mxArray* getField(const mxArray *pa, const char *field_name)
    {
        mxArray *field = mxGetField(pa, 0, field_name);
        if (! field)
            throw std::invalid_argument(std::string("Field ")+field_name+" not found in structure!");
        return field;
    }

    /*!
     * \brief Get char field
     *
     * \param pa MATLAB struct array
     * \param field_name field name
     * \return field string value
     */
    std::string getCharField(const mxArray *pa, const char *field_name)
    {
        mxArray *field = getField(pa, field_name);
        if (! mxIsChar(field))
            mexErrMsgIdAndTxt("psi:ptycho:arg:field:type:notchar", "Field %s not a char array!", field_name);
        if (mxGetNumberOfDimensions(field)!=2 || mxGetM(field)!=1)
            mexErrMsgIdAndTxt("psi:ptycho:arg:field:dim", "Field %s is not a single char array!", field_name);
        std::size_t nchars = mxGetN(field);
        char buf[nchars+1];
        if (mxGetString(field, buf, nchars+1))
            mexErrMsgIdAndTxt("psi:ptycho:string:extract", "Unable to extract characters of field %s!", field_name);
        return buf;
    }

    /*!
     * \brief Assign values to vector
     *
     * \param a Vector of length at least nelements
     * \param b Pointer to at least nelements values
     * \param nelements number of elements to assign
     * \tparam A Vector element type
     * \tparam B Value type
     */
    template <typename A, typename B>
    void assign(std::vector<A> &a, const B *b, std::size_t nelements)
    {
        if (b == nullptr)
            mexErrMsgIdAndTxt("psi:ptycho:arg:field:retrieval", "Unable to extract values from numeric vector field!");
        for (std::size_t i=0; i<nelements; i++)
            a[i] = b[i];
    }

    /*!
     * \brief Get vector field with long integer values
     *
     * \param pa MATLAB struct array
     * \param field_name field name
     * \param nelements expected number of elements (0 means undefined and is the default)
     * \return field string value
     * \tparam T Vector element type
     */
    template <typename T>
    std::vector<T> getVectorField(const mxArray *pa, const char *field_name, unsigned int nelements=0)
    {
        mxArray *field = getField(pa, field_name);
        if (! mxIsNumeric(field))
            mexErrMsgIdAndTxt("psi:ptycho:arg:field:type:notnumeric", "Field %s not a numeric array!", field_name);
        if (mxGetNumberOfDimensions(field)!=2 || mxGetM(field)!=1)
            mexErrMsgIdAndTxt("psi:ptycho:arg:field:dim", "Field %s is not a single numeric array!", field_name);
        std::size_t N = mxGetN(field);
        if (nelements && N != nelements)
            mexErrMsgIdAndTxt("psi:ptycho:arg:field:length", "Field %s doesn't have %u elements!", field_name, nelements);
        mxClassID mx_class = mxGetClassID(field);
        std::vector<T> buf(N);
        switch (mx_class) {
            case mxINT8_CLASS:
                assign(buf, mxGetInt8s(field), N);
                break;
            case mxUINT8_CLASS:
                assign(buf, mxGetUint8s(field), N);
                break;
            case mxINT16_CLASS:
                assign(buf, mxGetInt16s(field), N);
                break;
            case mxUINT16_CLASS:
                assign(buf, mxGetUint16s(field), N);
                break;
            case mxINT32_CLASS:
                assign(buf, mxGetInt32s(field), N);
                break;
            case mxUINT32_CLASS:
                assign(buf, mxGetUint32s(field), N);
                break;
            case mxINT64_CLASS:
                assign(buf, mxGetInt64s(field), N);
                break;
            case mxUINT64_CLASS:
                assign(buf, mxGetUint64s(field), N);
                break;
            case mxSINGLE_CLASS:
                assign(buf, mxGetSingles(field), N);
                break;
            case mxDOUBLE_CLASS:
                assign(buf, mxGetDoubles(field), N);
                break;
            default:
                std::string class_name(mxGetClassName(field));
                mexErrMsgIdAndTxt("psi:ptycho:arg:field:class", "Field %s doesn't contain integer elements (it has class %s)!", field_name, class_name.c_str());
        }
        return buf;
    }

    std::vector<std::string> getStrings(const mxArray *pa, const char *field_name, unsigned int nelements=0)
    {
        mxArray *field = getField(pa, field_name);
        if (mxGetM(field) != 1)
            mexErrMsgIdAndTxt("psi:ptycho:arg:wrongDimensions:notOne", "Argument %s must be a 1xN one-dimensional character or cell array!", field_name);
        std::unique_ptr<mx_trait<char>> ct(get_char_trait(field));
        if (ct.get()) { // char array
            return std::vector<std::string>({ ct->get_string(field) });
        } else {        // cell array
            if (! mxIsCell(field))
                mexErrMsgIdAndTxt("psi:ptycho:arg:illegal:type", "Argument %s must be either character or cell array!", field_name);
            std::vector<std::string> result_vec;
            for (std::size_t j=0; j<mxGetN(field); j+=1) {
                mxArray *pa = mxGetCell(field, j);
                if (mxGetM(pa) != 1)
                    mexErrMsgIdAndTxt("psi:ptycho:arg:wrongDimensions:notOne", "Cell fields in argument %s must be one-dimensional character arrays!", field_name);
                ct.reset(get_char_trait(pa));
                if (! ct.get())
                    mexErrMsgIdAndTxt("psi:ptycho:arg:wrongType:notChar", "Cell field in argument %s must be of type char!", field_name);
                result_vec.push_back(ct->get_string(pa));
                if (! result_vec.back().size())
                    mexErrMsgIdAndTxt("psi:ptycho:arg:illegal:empty", "Empty string in %s argument!", field_name);
            }
            if (! result_vec.size())
                mexErrMsgIdAndTxt("psi:ptycho:arg:illegal:empty", "Empty %s argument!", field_name);
            return result_vec;
        }
    }

} // namespace

/*!
 * \brief MEX function readMeasurementData
 *
 * Input arguments:
 * 0 : csax p structure
 *
 * Output arguments:
 * 0 : Measurement data (3d or 4d array [cols, rows, n_bursts, n_positions], n_burst is optional and only used if the HHDF5 dataset is 4d)
 *
 * \param nlhs Number of left hand side (result) arguments
 * \param plhs Pointer array to left hand side arguments
 * \param nrhs Number of right hand side (input) arguments
 * \param prhs Pointer array to right hand side arguments
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    constexpr unsigned int ninputs = 1;
    constexpr unsigned int noutputs = 1;

    DEBUG_INIT;
    if (nrhs != ninputs)
        mexErrMsgIdAndTxt("psi:ptycho:arg:in:wrongNumber", "Wrong number of arguments, use: array = readMeasurementData(argument structure)");
    if (nlhs != noutputs)
        mexErrMsgIdAndTxt("psi:ptycho:arg:out:wrongNumber", "Wrong number of output arguments, use: array = readObjectData(argument structure)");

    if (! mxIsStruct(prhs[0]))
        mexErrMsgIdAndTxt("psi:ptycho:arg:wrongType:notStruct", "Argument must be of type struct!");
    if (mxGetNumberOfDimensions(prhs[0])!=2 || mxGetM(prhs[0])!=1 || mxGetN(prhs[0])!=1)
        mexErrMsgIdAndTxt("psi:ptycho:arg:wrongDimensions:notOne", "Argument must be a 1x1 struct array!");

    bool verbose = env::is_enabled("PTYCHO_READ_VERBOSE");

    precision::type prec = precision::type::Single;
    try {
        prec = precision::from_str(getCharField(prhs[0], p::precision));
    } catch (std::exception &ex) {
        mexErrMsgIdAndTxt("psi:ptycho:arg:precision", "Precision field error: %s", ex.what());
    }
    try {
        std::string format(getCharField(prhs[0], p::format));
        std::vector<long> image_size;
        try {
            image_size = getVectorField<long>(prhs[0], p::image_size, 2);
        } catch (std::invalid_argument &ex) {}
        std::vector<long> roi_center;
        try {
            roi_center = getVectorField<long>(prhs[0], p::roi_center, 2);
        } catch (std::invalid_argument &ex) {}
        std::vector<long> nthreads(getVectorField<long>(prhs[0], p::nthreads, 1));
        std::vector<std::string> data_path(getStrings(prhs[0], p::data_path));
        // std::string scan_string_format(getCharField(prhs[0], p::scan_string_format));
        // std::string data_prefix(getCharField(prhs[0], p::data_prefix));
        // std::vector<long> scan_number(getVectorField<long>(prhs[0], p::scan_number));

        // Format dependent arguments
        std::vector<std::string> data_location;
        if (format == "h5")
            data_location = getStrings(prhs[0], p::data_location);

        if (verbose) {
            mexPrintf("extension: %s\n", format.c_str());
            if (data_path.size() == 1) {
                mexPrintf("data_path: %s\n", data_path[0].c_str());
            } else {
                mexPrintf("data_path first: %s\n", data_path[0].c_str());
                mexPrintf("data_path last : %s\n", data_path.back().c_str());
            }
            mexPrintf("roi_center: ");
            if (roi_center.empty()) {
                mexPrintf("from data\n");
            } else {
                mexPrintf("%ld %ld\n", roi_center[0], roi_center[1]);
            }
            mexPrintf("image_size: ");
            if (image_size.empty()) {
                mexPrintf("from data\n");
            } else {
                mexPrintf("%ld %ld\n", image_size[0], image_size[1]);
            }
            mexPrintf("precision: %s\n", precision::to_str(prec));
            //     mexPrintf("data_prefix: %s\n", data_prefix.c_str());
            //     mexPrintf("scan_number: [ ");
            //     for (long l : scan_number)
            //         mexPrintf("%ld ", l);
            //     mexPrintf("]\n");
            const char *debuglog = std::getenv("PTYCHO_READ_DEBUG");
            mexPrintf("debuglog: %s - ", (debuglog ? debuglog : "not set"));
            DEBUG {
                mexPrintf("ok\n");
                OUT << "----------- read measurement data ----------------\n";
            } else {
                mexPrintf("not functional\n");
            }

            // Format dependend verbose output
            if (format == "h5") {
                if (data_location.size() == 1) {
                    mexPrintf("data_location: %s\n", data_location[0].c_str());
                } else {
                    mexPrintf("data_location first: %s\n", data_location[0].c_str());
                    mexPrintf("data_location last : %s\n", data_location.back().c_str());
                }
            }
        }

        if (format == "h5") {
            DEBUG {
                OUT << "Reading Eiger data..." << std::endl;
            }
            plhs[0] = data_prep::read_eiger_data(nthreads[0], data_path, data_location, image_size, roi_center, prec);
        } else if (format == "cbf") {
            DEBUG {
                OUT << "Reading Pilatus data..." << std::endl;
            }
            plhs[0] = data_prep::read_data_threaded("pilatus", nthreads[0], data_path, image_size, roi_center, prec);
        } else if (format == "tiff") {
            DEBUG {
                OUT << "Reading Moench data..." << std::endl;
            }
            plhs[0] = data_prep::read_data_threaded("moench", nthreads[0], data_path, image_size, roi_center, prec);
        } else {
            mexErrMsgIdAndTxt("psi:ptycho:format:unknown", "Unknown extension: %s", format.c_str());
        }
    } catch (std::exception &ex) {
        mexErrMsgIdAndTxt("psi:ptycho:struct:arg", "%s", ex.what());
    }
    return;
}
