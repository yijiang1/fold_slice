/*!
 * \file
 * Helper definitions for HDF5 handling
 */

#include "hdf5_helper.h"

namespace {

    bool initialized = false;   //!< Datastructures initialized ?

} // namespace

namespace hdf5 {

    template<> hid_t type_trait<double>::type = H5T_NATIVE_DOUBLE;                  //!< HDF5 double
    template<> hid_t type_trait<float>::type = H5T_NATIVE_FLOAT;                    //!< HDF5 float
    template<> hid_t type_trait<uint64_t>::type = H5T_NATIVE_UINT64;                //!< HDF5 uint64_t
    template<> hid_t type_trait<unsigned int>::type = H5T_NATIVE_UINT;              //!< HDF5 unsigned int
    template<> hid_t type_trait<std::complex<float>>::type = type_trait<std::complex<float>>::invalid;      //!< HDF5 complex<float>
    template<> hid_t type_trait<std::complex<double>>::type = type_trait<std::complex<double>>::invalid;    //!< HDF5 complex<double>

    void initialize()
    {
        if (! initialized) {
            type_trait<std::complex<float>>::type = create_complex_type<double>();
            type_trait<std::complex<double>>::type = create_complex_type<float>();
        }
    }

    void cleanup()
    {
        if (initialized) {
            h5type_trait::close(type_trait<std::complex<float>>::type);
            h5type_trait::close(type_trait<std::complex<double>>::type);
        }
    }

}
