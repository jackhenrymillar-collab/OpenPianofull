/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

// Offline rendering of fixed test scenarios, used by the render tool and the regression test.

#include <cstdint>
#include <string>
#include <vector>

namespace openpiano::harness
{

struct Event
{
    enum Type { NoteOn, NoteOff, PedalDown, PedalUp, SilentKeyDown };

    int64_t sample;
    Type type;
    int note = 0;
    int velocity = 0;
};

struct Scenario
{
    std::string name;
    std::string description;
    double seconds;
    std::vector<Event> events; // Sorted by sample
    bool benchmarkOnly = false; // Excluded from the bit-exact regression test (slow)
    bool fullRangeOnly = false; // Uses keys the original 52-key OpenPiano doesn't have
};

// The standard set of scenarios, with event times computed for the given sample rate.
std::vector<Scenario> standardScenarios (int sampleRate);

// Renders through the restructured port of the original model (LegacyPiano), mono, splitting
// blocks at each event the way the plugin does.
std::vector<float> renderLegacy (const Scenario& scenario, int sampleRate, unsigned workerThreads, int blockSize);

// Renders through the new 88-key piano (PianoInstrument). Stereo, interleaved L R L R.
std::vector<float> renderPiano (const Scenario& scenario, int sampleRate, int blockSize);

// Renders through the original, unmodified OpenPiano code (tools/upstream), applying each
// event right before the sample it falls on.
std::vector<float> renderUpstream (const Scenario& scenario, int sampleRate);

bool writeWav (const std::string& path, const std::vector<float>& samples, int sampleRate, int channels = 1);

} // namespace openpiano::harness
