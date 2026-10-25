/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include "openpiano/Instrument.h"
#include "openpiano/ParallelFor.h"

#include <array>
#include <atomic>
#include <memory>
#include <vector>

namespace openpiano
{

namespace piano { class StringModel; }

// Player-facing controls for the piano. The defaults reproduce upstream OpenPiano exactly.
struct PianoControls
{
    float outputGain = 1.0f;      // Linear gain on top of the engine's internal level
    float dynamicsCurve = 1.0f;   // Velocity -> hammer speed exponent (1 = linear, >1 = softer touch)
    float hammerHardness = 1.0f;  // Scales felt stiffness (brighter above 1)
    float damperStrength = 1.0f;  // Scales how quickly released notes are damped
};

class PianoInstrument final : public Instrument
{
public:
    // numWorkerThreads: extra threads used to compute strings in parallel (0 = audio thread only).
    explicit PianoInstrument (unsigned numWorkerThreads = defaultWorkerThreads());
    ~PianoInstrument() override;

    static unsigned defaultWorkerThreads();

    const char* name() const override { return "Piano"; }

    void prepare (double sampleRate, int maxBlockSize) override;
    void reset() override;

    void noteOn (int midiNote, int velocity) override;
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
