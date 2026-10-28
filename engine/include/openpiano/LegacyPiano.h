/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include "openpiano/Instrument.h"
#include "openpiano/ParallelFor.h"
#include "openpiano/PianoControls.h"

#include <array>
#include <atomic>
#include <memory>
#include <vector>

namespace openpiano
{

namespace piano { class StringModel; }


// The original OpenPiano finite-difference piano (52 keys, one string per key), kept so the
// regression test can prove the port matches upstream and for before/after comparisons.
class LegacyPiano final : public Instrument
{
public:
    // numWorkerThreads: extra threads used to compute strings in parallel (0 = audio thread only).
    explicit LegacyPiano (unsigned numWorkerThreads = defaultWorkerThreads());
    ~LegacyPiano() override;

    static unsigned defaultWorkerThreads();

    const char* name() const override { return "Piano"; }

    void prepare (double sampleRate, int maxBlockSize) override;
    void reset() override;

    void noteOn (int midiNote, int velocity) override;
    using Instrument::noteOff;
    void noteOff (int midiNote) override;
    void setSustain (bool pedalDown) override;

    void process (float* const* outputs, int numChannels, int numSamples) override;

    int lowestNote() const override;
    int highestNote() const override;
    int activeVoiceCount() const override { return activeVoices.load (std::memory_order_relaxed); }

    // Audio thread only. Takes effect immediately.
    void setControls (const PianoControls& controls);

    // Hammer speed [m/s] for a MIDI velocity under the current dynamics curve.
    double hammerVelocityFor (int velocity) const;

    // Level that upstream OpenPiano applied in its plugin.
    static constexpr float kEngineGain = 150.0f;

private:
    void renderChunk (float* out, int numSamples);
    piano::StringModel* stringFor (int midiNote);

    std::unique_ptr<ParallelFor> pool;
    std::vector<std::unique_ptr<piano::StringModel>> strings; // index 0 = lowestNote()
    std::vector<std::vector<double>> stringBuffers;
    int chunkSize = 256;
    double sampleRate = 0;

    PianoControls controls;
    std::array<bool, 128> keyHeld {};
    bool sustainDown = false;

    std::atomic<int> activeVoices { 0 };
};

} // namespace openpiano
