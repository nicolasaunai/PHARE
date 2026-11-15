#ifndef PHARE_CORE_UTILITIES_PARALLEL_HPP
#define PHARE_CORE_UTILITIES_PARALLEL_HPP

#include <string>
#include <cstddef>
#include <cstdlib>
#include <algorithm>
#include <exception>

#if defined(_OPENMP)
#include <omp.h>
#endif

// Threads of one MPI process. Built without OpenMP everything here is serial.

namespace PHARE::core
{
namespace detail
{
    inline std::size_t& nbrThreads_()
    {
        static std::size_t n = 1;
        return n;
    }
} // namespace detail


// threads parallel_for may use on this process, 1 unless set
inline std::size_t nbrThreads()
{
    return detail::nbrThreads_();
}

inline void setNbrThreads([[maybe_unused]] std::size_t const n)
{
#if defined(_OPENMP)
    detail::nbrThreads_() = std::max<std::size_t>(n, 1);
#endif
}

// PHARE_THREADS, 1 if unset: an MPI run that does not ask for threads never starts any
inline void setNbrThreadsFromEnv()
{
    auto const* const env = std::getenv("PHARE_THREADS");
    setNbrThreads(env ? std::stoul(env) : 1);
}

// index of the calling thread among the threads of the enclosing parallel_for, 0 outside
inline std::size_t threadIndex()
{
#if defined(_OPENMP)
    return static_cast<std::size_t>(omp_get_thread_num());
#else
    return 0;
#endif
}

// threads a parallel_for over n indices runs on, threadIndex() is below this in fn
inline std::size_t nbrThreadsFor(std::size_t const n)
{
    return std::max<std::size_t>(std::min(n, nbrThreads()), 1);
}


// fn(i) has run for every i in [0, n) on return. Indices are handed one at a time to idle
// threads, so the order is not defined: fn(i) and fn(j) must not write the same data.
// From inside a parallel_for it runs on the calling thread, which keeps its threadIndex().
// The first exception thrown by fn is rethrown here, once all indices are done.
template<typename Fn>
void parallel_for(std::size_t const n, Fn&& fn)
{
#if defined(_OPENMP)
    auto const nthreads = nbrThreadsFor(n);
    if (nthreads > 1 and !omp_in_parallel())
    {
        std::exception_ptr error;

#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads)
        for (std::size_t i = 0; i < n; ++i)
        {
            try
            {
                fn(i);
            }
            catch (...)
            {
#pragma omp critical(phare_parallel_for_error)
                if (!error)
                    error = std::current_exception();
            }
        }

        if (error)
            std::rethrow_exception(error);
        return;
    }
#endif

    for (std::size_t i = 0; i < n; ++i)
        fn(i);
}

} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_PARALLEL_HPP
