/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include <atomic>
#include <cstdint>
#include <thread>
#include <vector>

namespace openpiano
{

// A small worker pool for the audio thread.
//
// run() splits `count` independent jobs across the workers and the calling thread, and
// returns when all of them are finished. The caller always takes part, so a job list
// completes even if no worker wakes up in time. run() does not allocate or lock.
class ParallelFor
{
public:
    // numWorkers extra threads are started; 0 means everything runs on the caller.
    explicit ParallelFor (unsigned numWorkers);
    ~ParallelFor();

    ParallelFor (const ParallelFor&) = delete;
    ParallelFor& operator= (const ParallelFor&) = delete;

    template <typename Fn>
    void run (int count, Fn& fn)
    {
        runImpl (count, [] (void* ctx, int i) { (*static_cast<Fn*> (ctx)) (i); }, &fn);
    }

    unsigned numWorkers() const { return (unsigned) threads.size(); }

private:
    using JobFn = void (*) (void*, int);

    void runImpl (int count, JobFn fn, void* ctx);
    void workerLoop();
    bool claimAndRun (uint32_t generation);

    static constexpr uint32_t closedIndex = 0xffffffffu;
    static uint64_t pack (uint32_t generation, uint32_t index) { return ((uint64_t) generation << 32) | index; }
    static uint32_t generationOf (uint64_t s) { return (uint32_t) (s >> 32); }
    static uint32_t indexOf (uint64_t s) { return (uint32_t) s; }

    // High 32 bits: job generation. Low 32 bits: next job index, or closedIndex between runs.
    std::atomic<uint64_t> state { pack (0, closedIndex) };
    std::atomic<int> remaining { 0 };
    std::atomic<bool> quit { false };

    // Written by the caller before a run is published, read by workers after claiming a job.
    std::atomic<int> jobCount { 0 };
    JobFn jobFn = nullptr;
    void* jobCtx = nullptr;

    std::vector<std::thread> threads;
};

} // namespace openpiano
