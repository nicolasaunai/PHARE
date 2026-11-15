#include "core/utilities/parallel.hpp"

#include "gtest/gtest.h"

#include <atomic>
#include <chrono>
#include <thread>
#include <algorithm>
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


TEST(ParallelFor, threadIndexIsBelowTheNumberOfThreads)
{
    Threads const threads{4};

    for (std::size_t const n : {std::size_t{2}, std::size_t{100}})
    {
        std::atomic<std::size_t> outOfRange{0};
        parallel_for(n, [&](std::size_t) {
            if (threadIndex() >= nbrThreads())
                ++outOfRange;
        });
        EXPECT_EQ(outOfRange, 0u);
    }
}


TEST(ParallelFor, nestedCallVisitsEveryIndexOnce)
{
    Threads const threads{4};
    std::size_t const nOuter = 3, nInner = 50;

    std::vector<int> visits(nOuter * nInner, 0);
    std::atomic<std::size_t> outOfRange{0};
    parallel_for(nOuter, [&](std::size_t const i) {
        parallel_for(nInner, [&](std::size_t const j) {
            ++visits[i * nInner + j];
            if (threadIndex() >= nbrThreads())
                ++outOfRange;
        });
    });

    for (auto const v : visits)
        EXPECT_EQ(v, 1);
    EXPECT_EQ(outOfRange, 0u);
}


TEST(ParallelFor, nestedCallIsSharedWithTheThreadsLeftWithoutAnOuterIndex)
{
    Threads const threads{4};

    // one outer index: three threads have nothing else to do than the inner indices
    std::vector<std::size_t> ranOn(32);
    parallel_for(1, [&](std::size_t) {
        parallel_for(ranOn.size(), [&](std::size_t const j) {
            ranOn[j] = threadIndex();
            volatile double x = 0;
            for (int k = 0; k < 2000000; ++k)
                x = x + k;
        });
    });

    std::vector<int> used(nbrThreads(), 0);
    for (auto const t : ranOn)
        used[t] = 1;
    EXPECT_GT(std::count(used.begin(), used.end(), 1), 1);
}


TEST(ParallelFor, manyNestedLoopsLoseNoIndex)
{
    Threads const threads{8};
    std::size_t const nOuter = 37, nInner = 23, nRepeats = 200;

    for (std::size_t r = 0; r < nRepeats; ++r)
    {
        std::atomic<std::size_t> sum{0};
        parallel_for(nOuter, [&](std::size_t const i) {
            parallel_for(nInner, [&](std::size_t const j) { sum += i * nInner + j + 1; });
            parallel_for(2, [&](std::size_t) { ++sum; });
        });
        auto const n = nOuter * nInner;
        ASSERT_EQ(sum, n * (n + 1) / 2 + 2 * nOuter);
    }
}


TEST(ParallelFor, workersThatWentToSleepComeBack)
{
    Threads const threads{4};

    std::this_thread::sleep_for(std::chrono::milliseconds{20}); // longer than they spin
    std::vector<int> visits(1000, 0);
    parallel_for(visits.size(), [&](std::size_t const i) { ++visits[i]; });
    for (auto const v : visits)
        EXPECT_EQ(v, 1);
}


TEST(ParallelFor, nestedCallRethrows)
{
    Threads const threads{4};

    EXPECT_THROW(parallel_for(2,
                              [&](std::size_t) {
                                  parallel_for(16, [&](std::size_t const j) {
                                      if (j == 5)
                                          throw std::runtime_error("inner 5");
                                  });
                              }),
                 std::runtime_error);
}


TEST(Chunks, coverTheRangeOnceInOrder)
{
    Threads const threads{4};

    EXPECT_EQ(nbrChunks(10, 100), 1u);
    EXPECT_EQ(nbrChunks(1000, 100), 10u);

    std::size_t const n = 1003, nChunks = nbrChunks(n, 100);
    std::size_t expectedFirst = 0;
    for (std::size_t c = 0; c < nChunks; ++c)
    {
        auto const [first, last] = chunk(n, nChunks, c);
        EXPECT_EQ(first, expectedFirst);
        EXPECT_GE(last - first, 100u);
        expectedFirst = last;
    }
    EXPECT_EQ(expectedFirst, n);
}


TEST(Chunks, areOneWithOneThread)
{
    EXPECT_EQ(nbrChunks(1000000, 100), 1u);
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
    EXPECT_EQ(nbrThreads(), 3u);

    unsetenv("PHARE_THREADS");
    setNbrThreadsFromEnv();
    EXPECT_EQ(nbrThreads(), 1u);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
