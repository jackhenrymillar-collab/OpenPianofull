/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

// Checks that the restructured engine renders exactly what upstream OpenPiano renders,
// sample for sample, across thread counts, block sizes and sample rates.
//
// Once the model is deliberately changed, this test will be switched to compare against
// stored renders of the new model instead.

#include "Harness.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <cstdio>
#include <cstring>

using namespace openpiano::harness;

int main (int argc, char** argv)
{
    // --quick runs a reduced matrix (used by CI on slow machines)
    const bool quick = argc > 1 && std::strcmp (argv[1], "--quick") == 0;

    struct Config { int rate; unsigned threads; int block; };
    const Config full[] = { { 48000, 0, 256 }, { 48000, 3, 256 }, { 48000, 2, 61 }, { 44100, 1, 512 } };
    const Config fast[] = { { 48000, 2, 256 } };

    int failures = 0;
    const Config* configs = quick ? fast : full;
    const size_t numConfigs = quick ? std::size (fast) : std::size (full);

    for (size_t c = 0; c < numConfigs; ++c)
    {
        const auto cfg = configs[c];
        for (auto& sc : standardScenarios (cfg.rate))
        {
            if (sc.benchmarkOnly)
                continue;

            const auto expected = renderUpstream (sc, cfg.rate);
            const auto actual = renderEngine (sc, cfg.rate, cfg.threads, cfg.block);

            size_t mismatches = 0, first = 0;
            float maxDiff = 0, peak = 0;
            for (size_t i = 0; i < expected.size(); ++i)
            {
                peak = std::max (peak, std::fabs (expected[i]));
                if (expected[i] != actual[i])
                {
                    if (mismatches++ == 0)
                        first = i;
                    maxDiff = std::max (maxDiff, std::fabs (expected[i] - actual[i]));
                }
            }

            const bool silent = peak == 0 && sc.name != "above_range_C6";
            const bool ok = mismatches == 0 && expected.size() == actual.size() && ! silent;
            std::printf ("%s  %5d Hz, %u workers, block %3d  %-18s",
                         ok ? "PASS" : "FAIL", cfg.rate, cfg.threads, cfg.block, sc.name.c_str());
            if (silent)
                std::printf ("  (upstream render is silent)");
            if (mismatches)
                std::printf ("  %zu samples differ, first at %zu, max diff %g", mismatches, first, maxDiff);
            std::printf ("\n");
            failures += ok ? 0 : 1;
        }
    }

    std::printf ("\n%s\n", failures == 0 ? "All renders match upstream bit for bit." : "MISMATCH");
    return failures == 0 ? 0 : 1;
}
