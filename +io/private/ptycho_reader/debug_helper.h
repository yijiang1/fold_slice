/*!
 * \file
 * Helper code for debug handling
 */

#ifndef DEBUG_HELPER
#define DEBUG_HELPER

/*!
 * \brief Debug functionality
 *
 * The environment variable PTYCHO_READ_DEBUG gives the name of the debug log file
 */
namespace debug {
    extern std::unique_ptr<std::ofstream> out; //!< pointer to debug file

    /*!
     * \brief Check debug output stream
     * \return true if stream is operational, otherwise false
     */
    inline bool check_out()
    {
        if (out.get())
            return out.get()->good();
        return false;
    }

    /*!
     * \brief Initialize debug output stream
     */
    inline void debug_init()
    {
        out.reset(nullptr);
        const char *fname = std::getenv("PTYCHO_READ_DEBUG");
        if (fname) {
            out.reset(new std::ofstream(fname, std::ios::app));
            if (! out.get()->good())
                out.reset(nullptr);
        }
    }
}

/*!
 * \brief Initialize debug stream
 */
#define DEBUG_INIT debug::debug_init()

/*!
 * \brief Start debug code block
 * Only access the debug stream inside such a block
 */
#define DEBUG if (debug::check_out())

/*!
 * \brief Access debug output stream
 * \return debug output stream reference
 */
#define OUT (*debug::out.get())

#endif
