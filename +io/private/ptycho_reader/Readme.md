This directory contains code for reading ptychographic reconstruction results and measurement data.

ptycho_read
===========

This is a MEX function for reading ptychographically constructed object datasets from HDF5 files.

The parameters to the function are
1. number of reader processes
2. precision of the returned array: either 'single' or 'double'
3. desired dimension of the returned objects
   This is an optional parameter. If absent, the dimension of the object in the first file is assumed to be the desired object dimension.
   Objects that have different dimension will be adapted (either some pixels at the borders left away, or zero pixels added)
4. path to the object dataset within the HDF5 file (e.g. '/reconstruction/object')
5. cell array with the file paths (e.g. {'the/path/to/example.h5'})
   If there is only one file, the parameter can also be a character array

From MATLAB call it in one of the following ways:

    [A, ind] = io.ptycho_read(1, 'single', [21,22], '/reconstruction/object', {'test1_20x20_c.h5', 'test2_20x20_c.h5'});
    A = io.ptycho_read(1, 'single', [21,22], '/reconstruction/object', {'test1_20x20_c.h5', 'test2_20x20_c.h5'});

The array *A* will contain the object datasets as a multidimensional array (Ncols, Nrows, Nmodes, Nslices, Nobjects). Nslices and Nmodes dimensions will be dropped if they are of size one. *ind* will contain a list of file indices that could not be read, if the first form is used. In the second form you'll get an error message if a file cannot be read. The function will print out some information about the read speed, if the environment variable *PTYCHO_READ_VERBOSE* is set to one of the values '1', 'yes', 'true'.

    setenv('PTYCHO_READ_VERBOSE', 'yes')

The above command accomplishes this inside MATLAB.

read_measurement
================

This is a MEX function for reading data from Eiger detector files in HDF5 format. The data must be located at /eh5/images within the file.

The parameters to the function are
1. parameter structure

The parameter structure contains some general fields:
* 'asize' with the result image dimensions (Ncols, Nrows) (optional, derived from data by default)
* 'ctr' with the regioin of interest center [col, row] relative to the data (optional, derived from data by default)
* 'precision' with the desired precision - 'single' or 'double'

For the **Eiger** detector, the parameter structure must contain the following fields:
* 'extension' with value 'h5'
* 'data_path' cell array (size 1 or N) with with the paths to the HDF5 files conaining the Eiger data
* 'data_location' cell array (size 1 or N) with the dataset location(s) within the file(s)
* 'nthreads' with the desired number of parallel read processes

Data is expected to be present in the following form:

    GROUP "/" {
       GROUP "eh5" {
          DATASET "images" {
             DATATYPE  H5T_STD_U32LE
             DATASPACE  SIMPLE { ( 464, 514, 1030 ) / ( 464, 514, 1030 ) }
          }
       }
    }

Here, there are 464 images with 514x1030 pixels. The file must be readable with the standard HDF5 library version 1.10.2

For the **Pilatus** detector
* 'extension' with value 'cbf'
* 'data_path' cell array with the paths of the CBF data files
* 'nthreads' with the desired number of parallel reader threads
Every data file must contain lines like

    conversions="x-CBF_BYTE_OFFSET"
    X-Binary-Size-Fastest-Dimension: 1475
    X-Binary-Size-Second-Dimension: 1679

where 1475x1679 are the image dimensions and the image data should be byte offset encoded.

For the **Moench** detector
* 'extension' with value 'tiff'
* 'data_path' cell array with the paths of the TIFF data files
* 'nthreads' with the desired number of parallel reader threads
The files are expected to contain raw single precision floating point data in one image plane, stored in IMAGELENGTH chunks of size IMAGEWIDTH. These tags must also give the image dimensions.

From MATLAB call it in the following way:

    arg = struct()
    arg.extension = 'h5'
    arg.precision = 'single'
    arg.nthreads = 1
    arg.data_path = { 'test.h5' }
    arg.data_location = { '/eiger/images' }
    A = io.read_measurement(arg);

The array *A* will contain the measurement datasets as a multidimensional array [Ncols, Nrows, Nimages]. The function will print out some information about the arguments, if the environment variable *PTYCHO_READ_VERBOSE* is set to one of the values '1', 'yes', 'true'.

    setenv('PTYCHO_READ_VERBOSE', 'yes')

The above command accomplishes this inside MATLAB.

Debug output
------------

If debug output about the inner workings of the functions is desired, the environment variable *PTYCHO_READ_DEBUG* should be set to the filename where debug output will end up.

    setenv('PTYCHO_READ_DEBUG', '/tmp/debug_output.txt')

The above command accomplishes this inside MATLAB.

Compilation
-----------

The code works only with MATLAB versions at or above 2018a, because the code requires the new interleaved complex array format. You need to load the matlab/2018a or a later MATLAB environment module in order to compile the code. Aditionally the gcc version and an apropriate HDF5 serial and TIFF environment module need to be loaded. This can be done by sourcing the *setup-environment.sh* script. After these steps module list shoud approximately like this:

    [stadler_h@ra-l-002 ~]$ module list
    Currently Loaded Modulefiles:
      1) gcc/6.3.0            2) hdf5_serial/1.8.18   3) matlab/2018a         4) tiff/4.0.9

Now you should be ready to compile the code:

    [stadler_h@ra-l-002 ptycho_reader]$ make
    mex -R2018a -c debug_helper.cc
    Building with 'g++'.
    MEX completed successfully.
    ...
    mex -R2018a readObjectData.cc read_object_data.o hdf5_helper.o debug_helper.o /opt/psi/Compiler/hdf5_serial/1.8.18/gcc/6.3.0/lib/libhdf5.a -lz -ldl -output ../../ptycho_read
    Building with 'g++'.
    MEX completed successfully.
    ...
    mex -R2018a readMeasurementData.cc read_eiger_data.o hdf5_helper.o read_data_threaded.o debug_helper.o /opt/psi/Compiler/hdf5_serial/1.8.18/gcc/6.3.0/lib/libhdf5.a -lz -ldl -output ../../read_measurement
    Building with 'g++'.
    MEX completed successfully.

The make process will procude the ptycho_read and the read_measurement MEX file in the io package: *io.ptycho_read* and *io.read_measurement*

Intermediate files of the make process can be cleaned using the command

    [stadler_h@ra-l-002 ptycho_reader]$ make clean
    rm -f *~ *.o

If you also want to delete the ptycho_read MEX file, use the command

    [stadler_h@ra-l-002 ptycho_reader]$ make proper
    rm -f *~ *.o
    rm -f ../../ptycho_read.mexa64 ../../read_measurement.mexa64

Code documentation in HTML format can be produced using doxygen (tested with version 1.8.13)

    [stadler_h@ra-l-002 ptycho_reader]$ make doc
    test -r doxygen.conf && rm -rf doc/html
    test -d doc || mkdir doc
    doxygen doxygen.conf
    ...

Hope everything works as expected!
