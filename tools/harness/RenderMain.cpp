/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

// openpiano-render: renders the standard test scenarios to WAV files and reports speed.
//
//   openpiano-render [--engine new|legacy|upstream|both] [--out DIR] [--rate HZ]
//                    [--threads N] [--block N] [--only TEXT] [--list]
//
//   new:      the 88-key modal piano (stereo)
//   legacy:   the restructured port of the original model (mono)
//   upstream: the original, unmodified OpenPiano code (mono)
//
// Files go to DIR/<engine>/<scenario>.wav (32-bit float). Compare two folders with
// tools/compare/compare.py.

#include "Harness.h"
#include "openpiano/LegacyPiano.h"
#include "piano/ScaleDesign.h"

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
    std::puts ("usage: openpiano-render [--engine new|legacy|upstream|both] [--out DIR] [--rate HZ]\n"
               "                        [--threads N] [--block N] [--only TEXT] [--list] [--scale]\n"
               "  --scale  print the 88-key string and hammer table as CSV");
}

int main (int argc, char** argv)
{
    std::string engine = "new", outDir = "renders", only;
    int rate = 48000, block = 256;
    unsigned threads = LegacyPiano::defaultWorkerThreads();
    bool list = false, scale = false;

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
        else if (arg == "--scale") scale = true;
        else { usage(); return arg == "--help" || arg == "-h" ? 0 : 2; }
    }

    if (engine != "new" && engine != "legacy" && engine != "upstream" && engine != "both") { usage(); return 2; }

    if (scale)
    {
        std::puts ("midi,strings,wound,length_m,core_mm,outer_mm,tension_N,B,f1_Hz,stretch_cents,hammer_g,Q0,p,strike_ratio,damper");
        for (const auto& k : piano::buildConcertGrandScale (0.0))
        {
            const auto& s = k.strings[0];
            std::printf ("%d,%d,%d,%.3f,%.3f,%.3f,%.0f,%.3g,%.3f,%.2f,%.2f,%.0f,%.3f,%.3f,%d\n", k.midiNote, k.numStrings,
                         s.wound ? 1 : 0, s.length, s.coreDiameter * 1e3, s.outerDiameter * 1e3, s.tension,
                         s.inharmonicity, k.targetFrequency, k.stretchCents, k.hammerMass * 1e3, k.feltStiffness,
                         k.feltExponent, k.strikeRatio, k.hasDamper ? 1 : 0);
        }
        return 0;
    }

    const auto scenarios = standardScenarios (rate);
    if (list)
    {
        for (auto& sc : scenarios)
            std::printf ("%-20s %s\n", sc.name.c_str(), sc.description.c_str());
        return 0;
    }

    std::printf ("Rate %d Hz, block %d, %u worker threads\n\n", rate, block, threads);
    std::printf ("%-9s %-20s %9s %10s %9s\n", "engine", "scenario", "peak dB", "render ms", "realtime");

    const char* engines[] = { "upstream", "legacy", "new" };
    for (const char* name : engines)
    {
        if (engine == "both" ? std::strcmp (name, "legacy") == 0 : engine != name)
            continue;
        const bool isNew = std::strcmp (name, "new") == 0;

        const auto dir = std::filesystem::path (outDir) / name;
        std::filesystem::create_directories (dir);

        for (auto& sc : scenarios)
        {
            if (! only.empty() && sc.name.find (only) == std::string::npos)
                continue;
            if (! isNew && sc.fullRangeOnly && sc.benchmarkOnly)
                continue;

            const auto start = std::chrono::steady_clock::now();
            const auto audio = isNew ? renderPiano (sc, rate, block)
                             : std::strcmp (name, "legacy") == 0 ? renderLegacy (sc, rate, threads, block)
                                                                 : renderUpstream (sc, rate);
            const int channels = isNew ? 2 : 1;
            const double ms = std::chrono::duration<double, std::milli> (std::chrono::steady_clock::now() - start).count();

            float peak = 0;
            for (float s : audio)
                peak = std::max (peak, std::fabs (s));
            const double peakDb = peak > 0 ? 20 * std::log10 (peak) : -INFINITY;

            // Realtime factor: how many times faster than playback (above 1 keeps up)
            const double rt = (sc.seconds * 1000.0) / ms;
            std::printf ("%-9s %-20s %9.1f %10.0f %8.1fx\n", name, sc.name.c_str(), peakDb, ms, rt);

            if (! writeWav ((dir / (sc.name + ".wav")).string(), audio, rate, channels))
            {
                std::fprintf (stderr, "could not write %s\n", (dir / (sc.name + ".wav")).string().c_str());
                return 1;
            }
        }
    }
    return 0;
}
