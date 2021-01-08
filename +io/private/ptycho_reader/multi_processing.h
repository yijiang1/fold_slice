#ifndef MULTI_PROCESSING_H
#define MULTI_PROCESSING_H

/*!
 * \file
 * Provide functionality for multi processing
 */

/*!
 * \brief Multiprocessing support functions and types
 */
namespace mp {

    /*!
     * \brief Buffer object managing destination and shared memory areas
     * \tparam e_type buffer element type
     */
    template <typename e_type>
    struct buf {
        const unsigned int msg_size;        //!< maximum message size per process
        e_type *destination;                //!< pointer to destination memory
        e_type *shared;                     //!< pointer to shared memory
        std::size_t length;                 //!< total data elements
        std::size_t chunk_size;             //!< data elements for one chunk
        int chunks_per_proc;                //!< number of chunks for a process to handle
        int first_full;                     //!< first process with non reduced number of chunks
        int num_procs;                      //!< number of read processes

        /*!
         * \brief Constructor
         *
         * \param dest Destination pointer
         * \param dims Result matrix dimensions
         * \param nprocs Number of read processes
         * \param message_size Size of message area per child process in bytes
         */
        buf(e_type *dest, const std::vector<mwSize> &dims, int nprocs, unsigned int message_size=64u)
            : msg_size(message_size ? message_size : 64u), destination(dest), shared(nullptr), length(0), chunk_size(1), chunks_per_proc(0), first_full(0), num_procs(nprocs)
        {
            assert(! dims.empty());
            assert(nprocs >= 1);
            int nchunks = dims[0];
            chunks_per_proc = (nchunks + nprocs - 1) / nprocs;
            first_full = (nprocs * chunks_per_proc) - nchunks;

            for (int i=1; i<dims.size(); i++)
                chunk_size *= dims[i];
            length = nchunks * chunk_size;

            DEBUG {
                OUT << nprocs << "np " << nchunks << "nc " << chunks_per_proc << "cpp " << first_full << "ff " << chunk_size << "cs " << length << 'l' << std::endl;
            }

            if (nprocs > 1) {
                void *sbuf = mmap(nullptr, (length - offset(1) * chunk_size) * sizeof(e_type) + (num_procs - 1) * msg_size, PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
                if (sbuf == MAP_FAILED)
                    throw std::runtime_error(std::string("mmap error - ") + strerror(errno));
                shared = static_cast<e_type *>(sbuf);
            }
        }

        /*!
         * \brief Destructor
         */
        ~buf() noexcept
        {
            if (shared)
                munmap(shared, (length - offset(1) * chunk_size) * sizeof(e_type) + (num_procs - 1) * msg_size);
        }

        /*!
         * \brief Get shareable buffer pointer
         *
         * \param proc process number
         * \return Pointer to buffer area accessible from parent process
         */
        e_type *get(int proc) noexcept
        {
            if (! proc) {
                DEBUG {
                    OUT << "buffer " << proc << ": destination+0 at address " << destination << std::endl;
                }
                return destination;
            }
            std::size_t buf_start = offset(proc, true) * chunk_size;
            DEBUG {
                OUT << "buffer " << proc << ": shared+" << buf_start << " at address " << &shared[buf_start] << std::endl;
            }
            return &shared[buf_start];
        }

        /*!
         * \brief Destination buffer
         *
         * \param proc process number
         * \return Pointer to destination buffer area for process proc
         */
        e_type *get_dest(int proc) noexcept
        {
            return destination + offset(proc) * chunk_size;
        }

        /*!
         * \brief Copy data in shared memory to destination
         */
        void copy() noexcept
        {
            if (shared) {
                DEBUG {
                    OUT << "Copy " << (length - offset(1) * chunk_size) << " elements to destination+" << (offset(1) * chunk_size) << std::endl;
                }
                std::copy(shared, shared + (length - offset(1) * chunk_size), destination + (offset(1) * chunk_size));
            }
        }

        /*!
         * \brief Object offset
         *
         * \param proc process number
         * \param shared offset into shared memory?
         * \return First chunk index for proc
         */
        int offset(int proc, bool shared=false) noexcept
        {
            int offset = proc * chunks_per_proc;
            if (proc < first_full)
                offset -= proc;
            else
                offset -= first_full;
            if (shared)
                offset -= chunks_per_proc - (first_full ? 1 : 0);
            return offset;
        }

        /*!
         * \brief Get pointer to message
         *
         * For child processes returns a pointer to the message area of size msg_size
         *
         * \param proc process number (must not be the parent process 0)
         * \return pointer to message area
         */
        char *get_msg(int proc)
        {
            if (proc)
                return reinterpret_cast<char *>(&shared[length - offset(1) * chunk_size]) + (proc - 1) * msg_size;
            throw std::logic_error("no message area for parent process");
            return nullptr;
        }

        /*!
         * \brief Record message
         *
         * Used to record a child process error message in the shared memory area
         *
         * \param msg message
         * \param proc process number (must not be the parent process 0)
         */
        void message(const char *msg, int proc)
        {
            std::strncpy(get_msg(proc), msg, msg_size-1);
        }
    }; // struct buf

    struct exception : public std::exception {
        bool converted = false;
        int process = 0;
        std::string error_message;

        inline exception ()
        {}

        inline exception (int proc, const std::string &msg)
            : process(proc), error_message(msg)
        {}

        inline const char* what() const noexcept override
        {
            try {
                if (! converted) {
                    std::string &msg = const_cast<std::string&>(error_message);
                    std::ostringstream oss;
                    oss << msg << " (proc " << process << ')';
                    msg = oss.str();
                    const_cast<bool&>(converted) = true;
                }
            } catch (...) {}
            return error_message.c_str();
        }
    }; // struct exception

    template<typename e_type>
    inline void run (int nprocs, buf<e_type> &buffer,
                     std::function<void(int, buf<e_type>&)> &&func)
    {
        if (nprocs <= 0)
            return;

        int proc = 0;             // process number starting with 0
        int pids[nprocs] = {0};   // first entry (parent) unused
        try {
            for (proc=nprocs-1; proc; --proc) {
                if ((pids[proc] = fork()) < 0) {
                    proc = 0;
                    throw exception(0, "process spawning failed");
                }
                if (! pids[proc])
                    break;
            }

            DEBUG_INIT;
            func(proc, buffer);

            if (proc) {
                DEBUG {
                    OUT << "Process " << proc << " finished" << std::endl;
                }
                std::_Exit((EXIT_SUCCESS));
            }
        } catch (std::exception &ex) {
            DEBUG {
                OUT << "Process " << proc << " exception: " << ex.what() << std::endl;
            }
            if (proc) {
                buffer.message(ex.what(), proc);
                std::_Exit((EXIT_FAILURE));
            }
            for (proc=nprocs-1; proc; --proc) {
                if (pids[proc] > 0) {
                    int status;
                    waitpid(pids[proc], &status, 0); // prevent zombies
                }
            }
            throw;
        }

        bool failed = false;
        exception ex;
        for (proc=nprocs-1; proc; --proc) {
            int status;
            bool with_error = false;
            if (waitpid(pids[proc], &status, 0) != pids[proc]) {
                with_error = true;
                ex.process = proc;
                ex.error_message = std::strerror(errno);
            } else if (! WIFEXITED(status)) {
                with_error = true;
                ex.process = proc;
                ex.error_message = "did not exit normally";
            } else if (WEXITSTATUS(status) != (EXIT_SUCCESS)) {
                with_error = true;
                ex.process = proc;
                ex.error_message = buffer.get_msg(proc);
            }
            DEBUG {
                OUT << "Process " << proc << ": finished " << (with_error ? "with error" : "successfully") << std::endl;
            }
            failed |= with_error;
        }
        if (failed)
            throw ex;

        buffer.copy();
    }

} // namespace mp

#endif // ifndef MULTI_PROCESSING_H
