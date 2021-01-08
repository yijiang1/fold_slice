/*!
 * \file
 * Helper code for HDF5 handling
 */

#include <complex>
#include <utility>
#include "hdf5.h"

/*!
 * \brief HDF5 helper classes and functions
 */
namespace hdf5 {

    /*!
     * \brief Type specific HDF5 features
     *
     * \tparam T native type
     */
    template <typename T>
    struct type_trait final {
        static constexpr hid_t invalid = -1;    //!< invalid type value
        static hid_t type;                      //!< HDF5 type handle for native type
    };

    /*!
     * \brief HDF5 object with id
     *
     * \tparam trait trait for HDF5 object type
     */
    template <typename trait>
    class obj final {
        obj(const obj &rhs) = delete;                           //!< no copy
        obj& operator= (const obj &rhs) = delete;               //!< no copy
        hid_t id;                                               //!< Wrapped object id
      public:
        static constexpr hid_t invalid = -1;                    //!< invalid id value
        /*!
         * \brief Constructor
         *
         * \param id_ set id to this
         */
        explicit obj(hid_t id_) noexcept : id(id_)
        {}

        /*!
         * \brief Destructor
         */
        ~obj() noexcept
        {
            if (valid())
                trait::close(id);
        }

        operator hid_t() noexcept { return id; }                //!< Conversion to hid_t
        bool valid() const noexcept { return id >= 0; }         //!< Valid? (true if id is valid)

        /*!
         * \brief Set new and return old id value
         *
         * \param val new id value
         * \return new id value
         */
        hid_t set(hid_t val) noexcept
        {
            if (valid())
                trait::close(id);
            id = val;
            return val;
        }

        /*!
         * \brief Grab id and invalidate this object
         *
         * \return id
         */
        hid_t grab() noexcept
        {
            hid_t val = id;
            id = invalid;
            return val;
        }
    };

    /*! \brief HDF5 File trait */
    struct file_trait final {
        static void close(hid_t id) noexcept { H5Fclose(id); }         //!< Close operation
    };

    /*! \brief HDF5 Attribute trait */
    struct attribute_trait final {
        static void close(hid_t id) noexcept { H5Aclose(id); }         //!< Close operation
    };

    /*! \brief HDF5 Group trait */
    struct group_trait final {
        static void close(hid_t id) noexcept { H5Gclose(id); }         //!< Close operation
    };

    /*! \brief HDF5 Dataset trait */
    struct dataset_trait final {
        static void close(hid_t id) noexcept { H5Dclose(id); }         //!< Close operation
    };

    /*! \brief HDF5 Dataspace trait */
    struct dataspace_trait final {
        static void close(hid_t id) noexcept { H5Sclose(id); }         //!< Close operation
    };

    /*! \brief HDF5 Type trait */
    struct h5type_trait final {
        static void close(hid_t id) noexcept { H5Tclose(id); }         //!< Close operation
    };

    /*! \brief HDF5 property trait */
    struct property_trait final {
        static void close(hid_t id) noexcept { H5Pclose(id); }         //!< Close operation
    };

    using file = obj<file_trait>;               //!< HDF5 file
    using attribute = obj<attribute_trait>;     //!< HDF5 attribute
    using group = obj<group_trait>;             //!< HDF5 group
    using dataset = obj<dataset_trait>;         //!< HDF5 dataset
    using dataspace = obj<dataspace_trait>;     //!< HDF5 dataspace
    using h5type = obj<h5type_trait>;           //!< HDF5 type
    using property = obj<property_trait>;       //!< HDF5 property

    /*!
     * \brief HDF5 exception
     */
    class exception final : public std::exception {
        std::string msg;    //!< error message

      public:
        /*!
         * \brief Constructor
         *
         * \param message error message
         */
        exception(std::string &&message)
            : msg(message)
        {}

        /*!
         * \brief Destructor
         */
        virtual ~exception() {}

        /*!
         * \brief Get error message
         *
         * \return error message
         */
        const char* what() const noexcept override
        {
            return msg.c_str();
        }
    };

    /*!
     * \brief Check hdf5 result
     *
     * \throw exception if result is not ok
     */
    inline void check_result(herr_t id, std::string &&message) {
        if (id < 0)
            throw exception(std::forward<decltype(message)>(message));
    }

    /*!
     * \brief Check hdf5 object
     *
     * \throw exception if object is not ok
     */
    inline void check_object(hid_t id, std::string &&message) {
        if (id < 0)
            throw exception(std::forward<decltype(message)>(message));
    }

    /*!
     * \brief Initialize data
     *
     * Must be called before using the complex data types
     */
    void initialize();

    /*!
     * \brief Cleanup data
     */
    void cleanup();

    /*!
     * \brief Convenience initializer class
     */
    struct initializer final {
        initializer() { initialize(); } //!< Constructor initializes internal hdf5 helper data structures
        ~initializer() { cleanup(); }   //!< Destructor cleans up internal hdf5 helper data structures
    };

    /*!
     * \brief Create complex compound type
     *
     * \return type id
     * \throw hdf5::exception if an error happens
     */
    template <typename T>
    hid_t create_complex_type()
    {
        hdf5::h5type type(H5Tcreate(H5T_COMPOUND, sizeof(std::complex<T>)));
        hdf5::check_result(type, "Unable to create complex compound type");
        hdf5::check_result(H5Tinsert(type, "r", 0, hdf5::type_trait<T>::type), "Unable to insert first complex compound member");
        hdf5::check_result(H5Tinsert(type, "i", sizeof(T), hdf5::type_trait<T>::type), "Unable to insert second complex compound member");
        return type.grab();
    }

} // namespace hdf5
