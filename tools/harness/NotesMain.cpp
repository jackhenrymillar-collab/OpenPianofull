/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

// openpiano-notes: renders single held notes (no pedal), one file per note and velocity,
// for comparing the model with recorded single-note samples (tools/calibrate).
//
//   openpiano-notes [--notes 21,33,45] [--vel 40,60,92,124] [--seconds S] [--out DIR]
//                   [--rate HZ] [--threads N]
//
// Each note is struck at time 0 and its key held until the end. Files are
// DIR/n<midi>_v<velocity>.wav (stereo, 32-bit float).

#include "Harness.h"
#include "openpiano/PianoInstrument.h"

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace openpiano;

static std::vector<int> parseList (const std::string& text)
{
    std::vector<int> out;
    std::stringstream ss (text);
    std::string item;
    while (std::getline (ss, item, ','))
        if (! item.empty())
            out.push_back (std::atoi (item.c_str()));
    return out;
}

int main (int argc, char** argv)
{
    std::vector<int> notes { 21, 33, 45, 57, 69, 81, 93, 105 };
    std::vector<int> velocities { 40, 60, 92, 124 };
    double seconds = 8.0;
    int rate = 48000;
    unsigned threads = std::max (1u, std::thread::hardware_concurrency());
    std::string outDir = "notes";

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        auto next = [&] () -> std::string
        {
            if (i + 1 >= argc) { std::fprintf (stderr, "missing value for %s\n", arg.c_str()); std::exit (2); }
            return argv[++i];
        };
        if (arg == "--notes") notes = parseList (next());
        else if (arg == "--vel") velocities = parseList (next());
        else if (arg == "--seconds") seconds = std::atof (next().c_str());
        else if (arg == "--out") outDir = next();
        else if (arg == "--rate") rate = std::atoi (next().c_str());
        else if (arg == "--threads") threads = (unsigned) std::max (1, std::atoi (next().c_str()));
        else
        {
            std::puts ("usage: openpiano-notes [--notes LIST] [--vel LIST] [--seconds S] [--out DIR] [--rate HZ] [--threads N]");
            return arg == "--help" || arg == "-h" ? 0 : 2;
        }
    }

    std::filesystem::create_directories (outDir);
    std::vector<std::pair<int, int>> jobs;
    for (int n : notes)
        for (int v : velocities)
            jobs.emplace_back (n, v);

    std::atomic<size_t> nextJob { 0 };
    std::atomic<bool> failed { false };
    auto worker = [&]
    {
        constexpr int block = 256;
        std::vector<float> left (block), right (block);
        for (size_t j; (j = nextJob++) < jobs.size();)
        {
            const auto [note, vel] = jobs[j];
            PianoInstrument piano;
            piano.prepare (rate, block);
            piano.noteOn (note, vel);
            const auto total = (size_t) (seconds * rate);
            std::vector<float> out (2 * total);
            for (size_t pos = 0; pos < total; pos += block)
            {
                const int n = (int) std::min<size_t> (block, total - pos);
                float* ch[] = { left.data(), right.data() };
                piano.process (ch, 2, n);
                for (int i = 0; i < n; ++i)
                {
                    out[2 * (pos + (size_t) i)] = left[(size_t) i];
                    out[2 * (pos + (size_t) i) + 1] = right[(size_t) i];
                }
            }
            const auto path = (std::filesystem::path (outDir) / ("n" + std::to_string (note) + "_v" + std::to_string (vel) + ".wav")).string();
            if (! harness::writeWav (path, out, rate, 2))
                failed = true;
        }
    };
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < std::min<unsigned> (threads, (unsigned) jobs.size()); ++t)
        pool.emplace_back (worker);
    for (auto& t : pool)
        t.join();
    if (failed)
    {
        std::fprintf (stderr, "could not write some files\n");
        return 1;
    }
    std::printf ("rendered %zu notes to %s\n", jobs.size(), outDir.c_str());
    return 0;
}
