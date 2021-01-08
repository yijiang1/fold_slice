ifndef MATLAB_HOME
  $(error "MATLAB_HOME is undefined, load matlab/2018a module!")
endif

ifndef HDF5_SERIAL_HOME
  $(error "HDF5_SERIAL_HOME is undefined, load gcc/6.3.0 and hdf5_serial/1.8.18 modules!")
endif

ifndef TIFF_INCLUDE_DIR
  ifndef TIFF_LIBRARY_DIR
    $(error "TIFF_INCLUDE_DIR or TIFF_LIBRARY_DIR undefined, load tiff/4.0.9 module!")
  endif
endif

MEXOPTIONS += -R2018a

IO_DIR := ../..
RDOBJ_TARGET_BASE := $(IO_DIR)/ptycho_read
RDMSR_TARGET_BASE := $(IO_DIR)/read_measurement
TARGET_SUFFIX := mexa64
RDOBJ_TARGET_FILE := $(RDOBJ_TARGET_BASE).$(TARGET_SUFFIX)
RDMSR_TARGET_FILE := $(RDMSR_TARGET_BASE).$(TARGET_SUFFIX)

all: $(RDOBJ_TARGET_FILE) $(RDMSR_TARGET_FILE)

ptycho_read: $(RDOBJ_TARGET_FILE)

read_measurement: $(RDMSR_TARGET_FILE)

debug_helper.o : debug_helper.cc debug_helper.h
	mex $(MEXOPTIONS) -c $<

hdf5_helper.o: hdf5_helper.cc hdf5_helper.h
	mex $(MEXOPTIONS) -c $<

read_object_data.o: read_object_data.cc read_object_data.h precision.h hdf5_helper.h debug_helper.h env_helper.h multi_processing.h
	mex $(MEXOPTIONS) -c $<

read_eiger_data.o: read_eiger_data.cc read_eiger_data.h precision.h hdf5_helper.h debug_helper.h env_helper.h multi_processing.h
	mex $(MEXOPTIONS) -c $<

read_data_threaded.o: read_data_threaded.cc read_data_threaded.h precision.h debug_helper.h env_helper.h ${TIFF_INCLUDE_DIR}/tiffio.h
	mex $(MEXOPTIONS) -I$(TIFF_INCLUDE_DIR) -c $<

# Put the resultinf mex file into +io/+HDF folder
$(RDOBJ_TARGET_FILE): readObjectData.cc debug_helper.o hdf5_helper.o read_object_data.o mex_helper.h env_helper.h precision.h
	mex $(MEXOPTIONS) $< read_object_data.o hdf5_helper.o debug_helper.o $(HDF5_SERIAL_HOME)/lib/libhdf5.a -lz -ldl -output $(RDOBJ_TARGET_BASE)

$(RDMSR_TARGET_FILE): readMeasurementData.cc read_eiger_data.o hdf5_helper.o read_data_threaded.o debug_helper.o env_helper.h mex_helper.h precision.h read_eiger_data.h read_data_threaded.h $(TIFF_LIBRARY_DIR)/libtiff.a
	mex $(MEXOPTIONS) $< read_eiger_data.o hdf5_helper.o read_data_threaded.o debug_helper.o $(HDF5_SERIAL_HOME)/lib/libhdf5.a $(TIFF_LIBRARY_DIR)/libtiff.a -lz -ldl -output $(RDMSR_TARGET_BASE)

doc:
	test -r doxygen.conf && rm -rf doc/html
	test -d doc || mkdir doc
	doxygen doxygen.conf

clean:
	rm -f *~ *.o

proper: clean
	rm -f $(RDOBJ_TARGET_FILE) $(RDMSR_TARGET_FILE)

.PHONY: all clean proper ptycho_read read_measurement doc
