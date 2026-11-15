#include "core/utilities/parallel.hpp"

#include "gtest/gtest.h"

#include <atomic>
#include <vector>
#include <cstddef>
#include <stdexcept>

using namespace PHARE::core;

namespace
{
// restores the serial default so tests do not depend on their order
struct Threads
{
    explicit Threads(std::size_t const n) { setNbrThreads(n); }
    ~Threads() { setNbrThreads(1); }
};
} // namespace


TEST(ParallelFor, isSerialByDefault)
{
    EXPECT_EQ(nbrThreads(), 1u);

    std::vector<std::size_t> order;
    parallel_for(5, [&](std::size_t const i) {
        EXPECT_EQ(threadIndex(), 0u);
        order.push_back(i);
    });
    EXPECT_EQ(order, (std::vector<std::size_t>{0, 1, 2, 3, 4}));
}


TEST(ParallelFor, visitsEveryIndexOnce)
{
    Threads const threads{4};
    std::size_t const n = 1000;

    std::vector<int> visits(n, 0);
    parallel_for(n, [&](std::size_t const i) { ++visits[i]; });

    for (auto const v : visits)
        EXPECT_EQ(v, 1);
}


TEST(ParallelFor, threadIndexIsBelowTheNumberOfThreadsUsed)
{
    Threads const threads{4};

    for (std::size_t const n : {std::size_t{2}, std::size_t{100}})
    {
        std::atomic<std::size_t> outOfRange{0};
        parallel_for(n, [&](std::size_t) {
            if (threadIndex() >= nbrThreadsFor(n))
                ++outOfRange;
        });
        EXPECT_EQ(outOfRange, 0u);
    }
    EXPECT_EQ(nbrThreadsFor(2), std::min<std::size_t>(2, nbrThreads()));
}


TEST(ParallelFor, nestedCallRunsOnTheCallingThread)
{
    Threads const threads{4};

    std::atomic<std::size_t> moved{0};
    parallel_for(8, [&](std::size_t) {
        auto const outer = threadIndex();
        parallel_for(8, [&](std::size_t) {
            if (threadIndex() != outer)
                ++moved;
        });
    });
    EXPECT_EQ(moved, 0u);
}


TEST(ParallelFor, rethrowsAnExceptionAfterAllIndicesRan)
{
    Threads const threads{4};
    std::size_t const n = 64;

    std::atomic<std::size_t> ran{0};
    EXPECT_THROW(parallel_for(n,
                              [&](std::size_t const i) {
                                  ++ran;
                                  if (i == 7)
                                      throw std::runtime_error("index 7");
                              }),
                 std::runtime_error);

    if (nbrThreads() > 1) // serial stops at the throw, as a plain loop does
    {
        EXPECT_EQ(ran, n);
    }
}


TEST(ParallelFor, emptyRangeCallsNothing)
{
    Threads const threads{4};
    parallel_for(0, [](std::size_t) { FAIL(); });
}


TEST(ParallelFor, threadsComeFromTheEnvironment)
{
    setenv("PHARE_THREADS", "3", 1);
    setNbrThreadsFromEnv();
#if defined(_OPENMP)
    EXPECT_EQ(nbrThreads(), 3u);
#else
    EXPECT_EQ(nbrThreads(), 1u);
#endif

    unsetenv("PHARE_THREADS");
    setNbrThreadsFromEnv();
    EXPECT_EQ(nbrThreads(), 1u);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
