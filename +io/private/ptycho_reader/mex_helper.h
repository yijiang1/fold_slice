#ifndef MEX_HELPER
#define MEX_HELPER

/*!
 * \file
 * Helper code for MEX data
 *
 * include \<cstdint\> before this one
 */

/*!
 * \brief Unique pointer for MATLAB allocated object
 *
 * \tparam T type of object
 */
template <typename T>
struct mx_deleter final {
    /*!
     * \brief deleter for pointer
     *
     * \param p pointer
     */
    void operator()(T *p) const noexcept {
        if (p)
            mxFree(p);
    }
};

/*!
 * \brief Unique pointer for MATLAB array
 */
template <>
struct mx_deleter<mxArray> final {
    /*!
     * \brief deleter for MATLAB array
     *
     * \param p MATLAB array
     */
    void operator()(mxArray *p) const noexcept {
        mxDestroyArray(p);
    }
};

/*!
 * \brief unique_ptr for MATLAB objects
 *
 * \tparam T object type
 */
template <typename T>
using mx_ptr = std::unique_ptr<T, mx_deleter<T>>;

/*!
 * \brief Type specific MATLAB features
 *
 * \tparam T native type
 */
template <typename T>
struct mx_trait final {
};

/*!
 * \brief Numeric real
 */
struct mx_numeric {
    /*!
     * \brief first int64_t
     * \param pa array
     * \param index array index
     * \return first int64_t in data
     */
    virtual int64_t get_int64(const mxArray *pa, std::size_t index=0) = 0;
};

/*!
 * \brief Double specific MATLAB features
 */
template <>
struct mx_trait<double> final : public mx_numeric {
    using type = double;                                    //!< base type
    using complex_type = mxComplexDouble;                   //!< complex type
    constexpr static mxClassID class_id = mxDOUBLE_CLASS;   //!< MATLAB class id

    /*!
     * \brief Get data pointer
     *
     * \param pa complex MATLAB array
     * \return pointer to complex array data
     */
    static complex_type* get_complex(mxArray *pa)
    {
        return mxGetComplexDoubles(pa);
    }
    /*!
     * \brief Get data pointer
     * \param pa MATLAB array
     * \return pointer to array data
     */
    static type* get(mxArray *pa)
    {
        return mxGetDoubles(pa);
    }
    /*!
     * \brief first int64_t
     * \param pa array
     * \param index array index
     * \return first int64_t in data
     */
    int64_t get_int64(const mxArray *pa, std::size_t index=0) override
    {
        return *(get(const_cast<mxArray*>(pa)) + index);
    }
};

/*!
 * \brief Float specific MATLAB features
 */
template <>
struct mx_trait<float> final : public mx_numeric {
    using type = float;                                     //!< base type
    using complex_type = mxComplexSingle;                   //!< complex type
    constexpr static mxClassID class_id = mxSINGLE_CLASS;   //!< MATLAB class id

    /*!
     * \brief Get data pointer
     *
     * \param pa complex MATLAB array
     * \return pointer to complex array data
     */
    static complex_type* get_complex(mxArray *pa)
    {
        return mxGetComplexSingles(pa);
    }
    /*!
     * \brief Get data pointer
     * \param pa MATLAB array
     * \return pointer to array data
     */
    static type* get(mxArray *pa)
    {
        return mxGetSingles(pa);
    }
    /*!
     * \brief first int64_t
     * \param pa array
     * \param index array index
     * \return first int64_t in data
     */
    int64_t get_int64(const mxArray *pa, std::size_t index=0) override
    {
        return *(get(const_cast<mxArray*>(pa)) + index);
    }
};

/*!
 * \brief int8 specific MATLAB features
 */
template <>
struct mx_trait<int8_t> final : public mx_numeric {
    using type = int8_t;                                    //!< base type
    constexpr static mxClassID class_id = mxINT8_CLASS;     //!< MATLAB class id
    /*!
     * \brief Get data pointer
     * \param pa MATLAB array
     * \return pointer to array data
     */
    static type* get(mxArray *pa)
    {
        return mxGetInt8s(pa);
    }
    /*!
     * \brief first int64_t
     * \param pa array
     * \param index array index
     * \return first int64_t in data
     */
    int64_t get_int64(const mxArray *pa, std::size_t index=0) override
    {
        return *(get(const_cast<mxArray*>(pa)) + index);
    }
};

/*!
 * \brief uint8 specific MATLAB features
 */
template <>
struct mx_trait<uint8_t> final : public mx_numeric {
    using type = uint8_t;                                   //!< base type
    constexpr static mxClassID class_id = mxUINT8_CLASS;    //!< MATLAB class id
    /*!
     * \brief Get data pointer
     * \param pa MATLAB array
     * \return pointer to array data
     */
    static type* get(mxArray *pa)
    {
        return mxGetUint8s(pa);
    }
    /*!
     * \brief first int64_t
     * \param pa array
     * \param index array index
     * \return first int64_t in data
     */
    int64_t get_int64(const mxArray *pa, std::size_t index=0) override
    {
        return *(get(const_cast<mxArray*>(pa)) + index);
    }
};

/*!
 * \brief int16 specific MATLAB features
 */
template <>
struct mx_trait<int16_t> final : public mx_numeric {
    using type = int16_t;                                   //!< base type
    constexpr static mxClassID class_id = mxINT16_CLASS;    //!< MATLAB class id
    /*!
     * \brief Get data pointer
     * \param pa MATLAB array
     * \return pointer to array data
     */
    static type* get(mxArray *pa)
    {
        return mxGetInt16s(pa);
    }
    /*!
     * \brief first int64_t
     * \param pa array
     * \param index array index
     * \return first int64_t in data
     */
    int64_t get_int64(const mxArray *pa, std::size_t index=0) override
    {
        return *(get(const_cast<mxArray*>(pa)) + index);
    }
};

/*!
 * \brief uint16 specific MATLAB features
 */
template <>
struct mx_trait<uint16_t> final : public mx_numeric {
    using type = uint16_t;                                  //!< base type
    constexpr static mxClassID class_id = mxUINT16_CLASS;   //!< MATLAB class id
    /*!
     * \brief Get data pointer
     * \param pa MATLAB array
     * \return pointer to array data
     */
    static type* get(mxArray *pa)
    {
        return mxGetUint16s(pa);
    }
    /*!
     * \brief first int64_t
     * \param pa array
     * \param index array index
     * \return first int64_t in data
     */
    int64_t get_int64(const mxArray *pa, std::size_t index=0) override
    {
        return *(get(const_cast<mxArray*>(pa)) + index);
    }
};

/*!
 * \brief int32 specific MATLAB features
 */
template <>
struct mx_trait<int32_t> final : public mx_numeric {
    using type = int32_t;                                   //!< base type
    constexpr static mxClassID class_id = mxINT32_CLASS;    //!< MATLAB class id
    /*!
     * \brief Get data pointer
     * \param pa MATLAB array
     * \return pointer to array data
     */
    static type* get(mxArray *pa)
    {
        return mxGetInt32s(pa);
    }
    /*!
     * \brief first int64_t
     * \param pa array
     * \param index array index
     * \return first int64_t in data
     */
    int64_t get_int64(const mxArray *pa, std::size_t index=0) override
    {
        return *(get(const_cast<mxArray*>(pa)) + index);
    }
};

/*!
 * \brief uint32 specific MATLAB features
 */
template <>
struct mx_trait<uint32_t> final : public mx_numeric {
    using type = uint32_t;                                  //!< base type
    constexpr static mxClassID class_id = mxUINT32_CLASS;   //!< MATLAB class id
    /*!
     * \brief Get data pointer
     * \param pa MATLAB array
     * \return pointer to array data
     */
    static type* get(mxArray *pa)
    {
        return mxGetUint32s(pa);
    }
    /*!
     * \brief first int64_t
     * \param pa array
     * \param index array index
     * \return first int64_t in data
     */
    int64_t get_int64(const mxArray *pa, std::size_t index=0) override
    {
        return *(get(const_cast<mxArray*>(pa)) + index);
    }
};

/*!
 * \brief int64 specific MATLAB features
 */
template <>
struct mx_trait<int64_t> final : public mx_numeric {
    using type = int64_t;                                   //!< base type
    constexpr static mxClassID class_id = mxINT64_CLASS;    //!< MATLAB class id
    /*!
     * \brief Get data pointer
     * \param pa MATLAB array
     * \return pointer to array data
     */
    static type* get(mxArray *pa)
    {
        return mxGetInt64s(pa);
    }
    /*!
     * \brief first int64_t
     * \param pa array
     * \param index array index
     * \return first int64_t in data
     */
    int64_t get_int64(const mxArray *pa, std::size_t index=0) override
    {
        return *(get(const_cast<mxArray*>(pa)) + index);
    }
};

/*!
 * \brief uint64 specific MATLAB features
 */
template <>
struct mx_trait<uint64_t> final : public mx_numeric {
    using type = uint64_t;                                  //!< base type
    constexpr static mxClassID class_id = mxUINT64_CLASS;   //!< MATLAB class id
    /*!
     * \brief Get data pointer
     * \param pa MATLAB array
     * \return pointer to array data
     */
    static type* get(mxArray *pa)
    {
        return mxGetUint64s(pa);
    }
    /*!
     * \brief first int64_t
     * \param pa array
     * \param index array index
     * \return first int64_t in data
     */
    int64_t get_int64(const mxArray *pa, std::size_t index=0) override
    {
        return *(get(const_cast<mxArray*>(pa)) + index);
    }
};

/*!
 * \brief get numeric trait for real numeric array
 * \param pa array
 * \return numeric trait or nullptr if not real numeric
 */
inline mx_numeric* get_numeric_trait(const mxArray *pa)
{
    if (mxIsComplex(pa))
        return nullptr;
    switch (mxGetClassID(pa)) {
        case mxDOUBLE_CLASS:
            return new mx_trait<double>;
        case mxSINGLE_CLASS:
            return new mx_trait<float>;
        case mxINT8_CLASS:
            return new mx_trait<int8_t>;
        case mxUINT8_CLASS:
            return new mx_trait<uint8_t>;
        case mxINT16_CLASS:
            return new mx_trait<int16_t>;
        case mxUINT16_CLASS:
            return new mx_trait<uint16_t>;
        case mxINT32_CLASS:
            return new mx_trait<int32_t>;
        case mxUINT32_CLASS:
            return new mx_trait<uint32_t>;
        case mxINT64_CLASS:
            return new mx_trait<int64_t>;
        case mxUINT64_CLASS:
            return new mx_trait<uint64_t>;
    }
    return nullptr;
}

template <>
struct mx_trait<char> final {
    constexpr static mxClassID class_id = mxCHAR_CLASS;     //!< MATLAB class id
    /*!
     * \brief Get string
     * \param pa array
     * \return string
     */
    static std::string get_string(const mxArray *pa)
    {
        const std::size_t nchars = mxGetN(pa);
        std::unique_ptr<char> buf(new char[nchars+1]);
        if (mxGetString(pa, buf.get(), nchars+1))
            return std::string();
        return std::string(buf.get(), nchars);
    }
};

/*!
 * \brief get char trait for class id
 * \param pa array
 * \return numeric trait or nullptr if not numeric
 */
inline mx_trait<char>* get_char_trait(const mxArray *pa)
{
    switch (mxGetClassID(pa)) {
        case mxCHAR_CLASS:
            return new mx_trait<char>;
    }
    return nullptr;
}

#endif // ifndef MEX_HELPER
