/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#include "openpiano/ParallelFor.h"

#include <chrono>

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__)
 #include <immintrin.h>
 #define OPENPIANO_CPU_RELAX() _mm_pause()
#elif defined(__aarch64__) || defined(__arm__)
 #define OPENPIANO_CPU_RELAX() __asm__ __volatile__ ("yield")
#else
 #define OPENPIANO_CPU_RELAX() ((void) 0)
#endif

namespace openpiano
{

ParallelFor::ParallelFor (unsigned numWorkers)
{
    threads.reserve (numWorkers);
    for (unsigned i = 0; i < numWorkers; ++i)
        threads.emplace_back ([this] { workerLoop(); });
}

ParallelFor::~ParallelFor()
{
    quit.store (true, std::memory_order_release);
    for (auto& t : threads)
        t.join();
}

bool ParallelFor::claimAndRun (uint32_t generation)
{
    bool didWork = false;
    uint64_t s = state.load (std::memory_order_acquire);

    for (;;)
    {
        if (generationOf (s) != generation || indexOf (s) == closedIndex)
            return didWork;

        const uint32_t index = indexOf (s);
        if ((int) index >= jobCount.load (std::memory_order_relaxed))
            return didWork;

        // The CAS only succeeds while this run is still open, which guarantees that
        // jobFn, jobCtx and jobCount all belong to this run.
        if (state.compare_exchange_weak (s, s + 1, std::memory_order_acq_rel, std::memory_order_acquire))
        {
            jobFn (jobCtx, (int) index);
            remaining.fetch_sub (1, std::memory_order_acq_rel);
            didWork = true;
            s = state.load (std::memory_order_acquire);
        }
    }
}

void ParallelFor::workerLoop()
{
    uint32_t lastGeneration = generationOf (state.load (std::memory_order_acquire));
    int idleSpins = 0;

    while (! quit.load (std::memory_order_acquire))
    {
        const uint64_t s = state.load (std::memory_order_acquire);
        const uint32_t g = generationOf (s);

        if (g != lastGeneration && indexOf (s) != closedIndex)
        {
            claimAndRun (g);
            lastGeneration = g;
            idleSpins = 0;
            continue;
        }

        // Spin briefly right after a block (the next one is usually a few ms away),
        // then back off to short sleeps so idle workers don't burn a core.
        if (++idleSpins < 4000)
            OPENPIANO_CPU_RELAX();
        else
            std::this_thread::sleep_for (std::chrono::microseconds (50));
    }
}

void ParallelFor::runImpl (int count, JobFn fn, void* ctx)
{
    if (count <= 0)
        return;

    if (threads.empty())
    {
        for (int i = 0; i < count; ++i)
            fn (ctx, i);
        return;
    }

    jobFn = fn;
    jobCtx = ctx;
    jobCount.store (count, std::memory_order_relaxed);
    remaining.store (count, std::memory_order_relaxed);

    const uint32_t generation = generationOf (state.load (std::memory_order_relaxed)) + 1;
    state.store (pack (generation, 0), std::memory_order_release);

    claimAndRun (generation);

    while (remaining.load (std::memory_order_acquire) > 0)
        OPENPIANO_CPU_RELAX();

    // Close the run before the next one rewrites jobFn/jobCtx/jobCount.
    state.store (pack (generation, closedIndex), std::memory_order_release);
}

} // namespace openpiano
