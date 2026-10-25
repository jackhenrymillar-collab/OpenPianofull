/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include "openpiano/Instrument.h"
#include "openpiano/PianoControls.h"

#include <array>
#include <atomic>
#include <memory>
#include <vector>

namespace openpiano
{

namespace piano { class ModalPiano; }

// The full 88-key piano: about 230 strings tuned as on a concert grand, hammers, a shared
// bridge and a soundboard, with sympathetic resonance. See piano/ModalPiano.h.
class PianoInstrument final : public Instrument
{
public:
    PianoInstrument();
    ~PianoInstrument() override;

    const char* name() const override { return "Piano"; }

    void prepare (double sampleRate, int maxBlockSize) override;
    void reset() override;

    void noteOn (int midiNote, int velocity) override;
    void noteOff (int midiNote) override;
    void setSustain (bool pedalDown) override;

    // Lifts a key's dampers without striking (a key pressed too gently to sound)
    void holdKeySilently (int midiNote);

    // Writes stereo to channels 0 and 1 (mono hosts get the sum on channel 0)
    void process (float* const* outputs, int numChannels, int numSamples) override;

    int lowestNote() const override;
    int highestNote() const override;
    int activeVoiceCount() const override { return activeStrings.load (std::memory_order_relaxed); }

    // Audio thread only. Takes effect immediately.
    void setControls (const PianoControls& controls);

    // Hammer speed [m/s] for a MIDI velocity under the current dynamics curve.
    double hammerVelocityFor (int velocity) const;

    // For tests and the render tool
    const piano::ModalPiano& model() const { return *piano; }

private:
    std::unique_ptr<piano::ModalPiano> piano;
    std::vector<float> scratch;
    PianoControls controls;
    std::array<bool, 128> keyHeld {};
    std::atomic<int> activeStrings { 0 };

    // Age changes the piano's build (tuning, strings, board), so it is applied by rebuilding
    // the model once the knob has rested for a moment
    double sampleRate = 48000;
    float builtAge = 0.0f;
    bool pedal = false;
    int ageSettleSamples = 0;
    void applyLiveControls();
};

} // namespace openpiano
