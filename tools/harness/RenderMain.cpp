/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

// openpiano-render: renders the standard test scenarios to WAV files and reports speed.
//
//   openpiano-render [--engine new|upstream|both] [--out DIR] [--rate HZ]
//                    [--threads N] [--block N] [--only TEXT] [--list]
//
// Files go to DIR/<engine>/<scenario>.wav (32-bit float, mono). Compare two folders with
// tools/compare/compare.py.

#include "Harness.h"
#include "openpiano/PianoInstrument.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <string>

using namespace openpiano;
using namespace openpiano::harness;

static void usage()
{
    std::puts ("usage: openpiano-render [--engine new|upstream|both] [--out DIR] [--rate HZ]\n"
               "                        [--threads N] [--block N] [--only TEXT] [--list]");
}

int main (int argc, char** argv)
{
    std::string engine = "new", outDir = "renders", only;
    int rate = 48000, block = 256;
    unsigned threads = PianoInstrument::defaultWorkerThreads();
    bool list = false;

    for (int i = 1; i < argc; ++i)
    {
        auto arg = std::string (argv[i]);
        auto next = [&] () -> std::string
        {
            if (i + 1 >= argc) { usage(); std::exit (2); }
            return argv[++i];
        };
        if (arg == "--engine") engine = next();
        else if (arg == "--out") outDir = next();
        else if (arg == "--rate") rate = std::atoi (next().c_str());
        else if (arg == "--threads") threads = (unsigned) std::atoi (next().c_str());
        else if (arg == "--block") block = std::atoi (next().c_str());
        else if (arg == "--only") only = next();
        else if (arg == "--list") list = true;
        else { usage(); return arg == "--help" || arg == "-h" ? 0 : 2; }
    }

    if (engine != "new" && engine != "upstream" && engine != "both") { usage(); return 2; }

    const auto scenarios = standardScenarios (rate);
    if (list)
    {
        for (auto& sc : scenarios)
            std::printf ("%-20s %s\n", sc.name.c_str(), sc.description.c_str());
        return 0;
    }

    std::printf ("Rate %d Hz, block %d, %u worker threads\n\n", rate, block, threads);
    std::printf ("%-9s %-20s %9s %10s %9s\n", "engine", "scenario", "peak dB", "render ms", "realtime");

    const char* engines[] = { "upstream", "new" };
    for (const char* name : engines)
    {
        if (engine != "both" && engine != name)
            continue;

        const auto dir = std::filesystem::path (outDir) / name;
        std::filesystem::create_directories (dir);

        for (auto& sc : scenarios)
        {
            if (! only.empty() && sc.name.find (only) == std::string::npos)
                continue;

            const auto start = std::chrono::steady_clock::now();
            const auto audio = std::strcmp (name, "new") == 0 ? renderEngine (sc, rate, threads, block)
                                                               : renderUpstream (sc, rate);
            const double ms = std::chrono::duration<double, std::milli> (std::chrono::steady_clock::now() - start).count();

            float peak = 0;
            for (float s : audio)
                peak = std::max (peak, std::fabs (s));
            const double peakDb = peak > 0 ? 20 * std::log10 (peak) : -INFINITY;

            // Realtime factor: how many times faster than playback (above 1 keeps up)
            const double rt = (sc.seconds * 1000.0) / ms;
            std::printf ("%-9s %-20s %9.1f %10.0f %8.1fx\n", name, sc.name.c_str(), peakDb, ms, rt);

            if (! writeWav ((dir / (sc.name + ".wav")).string(), audio, rate))
            {
                std::fprintf (stderr, "could not write %s\n", (dir / (sc.name + ".wav")).string().c_str());
                return 1;
            }
        }
    }
    return 0;
}
