/*!
 * \file
 * Read Object Data
 *
 * This file contains the main mex code for reading object data into MATLAB memory.
 */

#include <string>
#include <cstdlib>
#include <cstdint>
#include <memory>
#include <vector>
#include <stdexcept>
#include <fstream>
#include "precision.h"
#include "env_helper.h"
#include "mex.h"
#include "mex_helper.h"
#include "read_object_data.h"
#include "debug_helper.h"

#define USAGE "use: array = readObjectData(nProcesses, precision, [nrows, ncols], datasetPath, { filePaths, ... })"

/*!
 * \brief MEX function readObjectData
 *
 * Input arguments:
 * 0 : Number of read processes
 * 1 : 'double' or 'single' precision
 * 2 : Output image size [nrows, ncolumns] (optional, default value = size of first image)
 * 3 : Input object dataset path within the HDF5 file
 * 4 : Input HDF5 file paths
 *
 * Output arguments:
 * 0 : Object data (3d, 4d[with modes], or 5d[with slices]) [nobjects, nslices*, nmodes*, nrows, ncolumns]
 * 1 : (Optional) 1d array of file indices that could not be read
 *
 * \param nlhs Number of left hand side (result) arguments
 * \param plhs Pointer array to left hand side arguments
 * \param nrhs Number of right hand side (input) arguments
 * \param prhs Pointer array to right hand side arguments
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    constexpr int max_nprocs = 512;
    constexpr unsigned int min_inputs = 4;
    constexpr unsigned int max_inputs = 5;
    constexpr unsigned int min_outputs = 1;
    constexpr unsigned int max_outputs = 2;

    DEBUG_INIT;
    if ((nrhs < min_inputs) || (nrhs > max_inputs))
        mexErrMsgIdAndTxt("psi:ptycho:arg:in:wrongNumber", "Wrong number of arguments, " USAGE);
    if ((nlhs < min_outputs) || (nlhs > max_outputs))
        mexErrMsgIdAndTxt("psi:ptycho:arg:out:wrongNumber", "Wrong number of output arguments, " USAGE);

    unsigned int i = 0;
    int64_t nprocs = 1;
    { // 0
        if (! mxIsScalar(prhs[i]))
            mexErrMsgIdAndTxt("psi:ptycho:arg:wrongType:notScalar", "Argument %d must be scalar!", i+1);
        std::unique_ptr<mx_numeric> nt(get_numeric_trait(prhs[i]));
        if (! nt.get())
            mexErrMsgIdAndTxt("psi:ptycho:arg:wrongType:notScalar", "Argument %d must be numeric!", i+1);
        nprocs = nt->get_int64(prhs[i]);
        if (nprocs < 1)
            mexErrMsgIdAndTxt("psi:ptycho:arg:illegal:notPositive", "Number of processes must be positive!");
        if (nprocs > max_nprocs)
            mexErrMsgIdAndTxt("psi:ptycho:arg:illegal:tooMany", "Too many processes, max is %d!", max_nprocs);
        i += 1;
    }
    precision::type prec = precision::type::Single;
    { // 1
        if (mxGetM(prhs[i]) != 1)
            mexErrMsgIdAndTxt("psi:ptycho:arg:wrongDimensions:notOne", "Argument %d must be one-dimensional character array!", i+1);
        std::unique_ptr<mx_trait<char>> ct(get_char_trait(prhs[i]));
        if (! ct.get())
            mexErrMsgIdAndTxt("psi:ptycho:arg:wrongType:notChar", "Argument %d must be of type char!", i+1);
        std::string precString = ct->get_string(prhs[i]);
        if (precString == "double") {
            prec = precision::type::Double;
        } else if (precString != "single") {
            mexErrMsgIdAndTxt("psi:ptycho:arg:type", "Illegal type argument (%s): must be either 'single' or 'double'!", precString.c_str());
        }
        i += 1;
    }
    std::array<int64_t, 2> dims{0, 0};
    do { // 2 optionial
        std::unique_ptr<mx_numeric> nt(get_numeric_trait(prhs[i]));
        if (! nt.get())
            break;
        if ((mxGetN(prhs[i]) != 2) || (mxGetM(prhs[i]) != 1))
            mexErrMsgIdAndTxt("psi:ptycho:arg:wrongDimensions:oneByTwo", "Argument %d must be 2x1 numeric array!", i+1);
        dims[0] = nt->get_int64(prhs[i], 0);
        dims[1] = nt->get_int64(prhs[i], 1);
        if ((dims[0] < 1) || (dims[1] < 1))
            mexErrMsgIdAndTxt("psi:ptycho:arg:illegal:illegalSize", "Size argument must be at least 1 in every dimension!");
        i += 1;
    } while(false);
    std::string object_path;
    { // 3
        if (mxGetM(prhs[i]) != 1)
            mexErrMsgIdAndTxt("psi:ptycho:arg:wrongDimensions:notOne", "Argument %d must be one-dimensional character array!", i+1);
        std::unique_ptr<mx_trait<char>> ct(get_char_trait(prhs[i]));
        if (! ct.get())
            mexErrMsgIdAndTxt("psi:ptycho:arg:wrongType:notChar", "Argument %d must be of type char!", i+1);
        object_path = ct->get_string(prhs[i]);
        if (! object_path.size())
            mexErrMsgIdAndTxt("psi:ptycho:arg:illegal:empty", "Illegal object path argument!");
        i += 1;
    }
    if (i == nrhs)
        mexErrMsgIdAndTxt("psi:ptycho:arg:in:wrongNumber", "Missing arguments, " USAGE);
    std::vector<std::string> file_paths;
    { // 4
        if (mxGetM(prhs[i]) != 1)
            mexErrMsgIdAndTxt("psi:ptycho:arg:wrongDimensions:notOne", "Argument %d must be a 1xN one-dimensional character or cell array!", i+1);
        std::unique_ptr<mx_trait<char>> ct(get_char_trait(prhs[i]));
        if (ct.get()) { // char array
            file_paths.push_back(ct->get_string(prhs[i]));
        } else {        // cell array
            if (! mxIsCell(prhs[i]))
                mexErrMsgIdAndTxt("psi:ptycho:arg:illegal:type", "Argument %d must be either character or cell array!", i+1);
            for (std::size_t j=0; j<mxGetN(prhs[i]); j+=1) {
                mxArray *pa = mxGetCell(prhs[i], j);
                if (mxGetM(pa) != 1)
                    mexErrMsgIdAndTxt("psi:ptycho:arg:wrongDimensions:notOne", "Cell fields in argument %d must be one-dimensional character arrays!", i+1);
                ct.reset(get_char_trait(pa));
                if (! ct.get())
                    mexErrMsgIdAndTxt("psi:ptycho:arg:wrongType:notChar", "Cell field in argument %d must be of type char!", i+1);
                file_paths.push_back(ct->get_string(pa));
                if (! file_paths.back().size())
                    mexErrMsgIdAndTxt("psi:ptycho:arg:illegal:empty", "Illegal file path argument!");
            }
            if (! file_paths.size())
                mexErrMsgIdAndTxt("psi:ptycho:arg:illegal:empty", "Empty file path argument!");
        }
        i += 1;
    }
    if (i != nrhs)
        mexErrMsgIdAndTxt("psi:ptycho:arg:in:wrongNumber", "Superfluous arguments, " USAGE);

    bool verbose = env::is_enabled("PTYCHO_READ_VERBOSE");
    if (verbose) {
        mexPrintf("precision: %s\n", precision::to_str(prec));
        const char *debuglog = std::getenv("PTYCHO_READ_DEBUG");
        mexPrintf("debuglog: %s - ", (debuglog ? debuglog : "not set"));
        DEBUG {
            mexPrintf("ok\n");
            OUT << "----------- read object data ----------------\n";
        } else {
            mexPrintf("not functional\n");
        }
        mexPrintf("object_location: %s\n", object_path.c_str());
        if (file_paths.size() == 1) {
            mexPrintf("data_location: %s\n", file_paths[0].c_str());
        } else {
            mexPrintf("data_location first: %s\n", file_paths[0].c_str());
            mexPrintf("data_location last : %s\n", file_paths.back().c_str());
        }
    }

    data_prep::read_stat r_stat;
    plhs[0] = data_prep::read_object_data((int)nprocs, prec, dims, object_path, file_paths, (nlhs == 2) ? &plhs[1] : nullptr, verbose ? &r_stat : nullptr);
    if (verbose)
        mexPrintf("Read %fMB in %fs, bandwidth is %fMB/s\n", (r_stat.nbytes / 1000000.), r_stat.seconds, (r_stat.nbytes / (r_stat.seconds * 1000000.)));
    return;
}
