#ifndef PRECISION_H
#define PRECISION_H

/*!
 * \file
 * Precision handling support
 *
 * \pre needs standard exceptions and string
 */

/*!
 * \brief Precision related functions and types
 */
namespace precision {

    /*!
     * \brief Single or double precision
     */
    enum type : unsigned char {
        Single,    //!< single precision
        Double     //!< double precision
    };

    /*!
     * \brief Precision from string
     * \arg str precision string
     * \return precision
     */
    inline type from_str(const std::string &str)
    {
        type prec = Single;
        if (str == "double")
            prec = Double;
        else if (str != "single")
            throw std::invalid_argument("undefined precision string, must be 'single' or 'double'");
        return prec;
    }

    /*!
     * \brief Precision to string
     * \arg prec precision
     * \return string
     */
    inline const char* to_str(type prec)
    {
        switch (prec) {
            case Single:
                return "single";
            case Double:
                return "double";
            default:
                throw std::invalid_argument("undefined precision value");
        }
        return "undef";
    }

} // namespace precision

#endif
