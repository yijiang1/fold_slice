#ifndef ENV_HELPER_H
#define ENV_HELPER_H

/*!
 * \file
 * Helper functions for environment parsing
 */

namespace env {

    /*!
     * \brief Check if environment variable steered feature is enabled
     *
     * \param envvar environment variable name
     * \return true if feature is enabled
     */
    inline bool is_enabled (const std::string &envvar)
    {
        char *val = std::getenv(envvar.c_str());
        if (! val)
            return false;
        std::string value(val);
        if (value == "1")
            return true;
        if (value == "yes")
            return true;
        if (value == "true")
            return true;
        return false;
    }

} // namespace

#endif
