/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#include "openpiano/PianoInstrument.h"

#include "ModalPiano.h"

#include <algorithm>
#include <cmath>

namespace openpiano
{

namespace
{
    // Leaves normal playing untouched; rounds off only the rare peaks above -2 dBFS,
    // so a fistful of fortissimo notes with the pedal down never clips hard.
    inline float safetyClip (float x)
    {
        constexpr float knee = 0.8f, range = 0.2f;
        const float ax = std::fabs (x);
        if (ax <= knee)
            return x;
        return std::copysign (knee + range * std::tanh ((ax - knee) / range), x);
    }
}

PianoInstrument::PianoInstrument() : piano (std::make_unique<piano::ModalPiano>()) {}
PianoInstrument::~PianoInstrument() = default;

int PianoInstrument::lowestNote() const { return piano::kLowestKey; }
int PianoInstrument::highestNote() const { return piano::kHighestKey; }

void PianoInstrument::prepare (double sampleRate, int maxBlockSize)
{
    piano->prepare (sampleRate, {});
    piano->setHammerHardness (controls.hammerHardness);
    piano->setDamperStrength (controls.damperStrength);
    piano->setActionNoise (controls.actionNoise);
    scratch.assign ((size_t) std::max (1, maxBlockSize), 0.0f);
    keyHeld.fill (false);
    activeStrings = 0;
}

void PianoInstrument::reset()
{
    piano->reset();
    keyHeld.fill (false);
}

double PianoInstrument::hammerVelocityFor (int velocity) const
{
    // Hammer speed from 0.8 m/s (pianissimo) to 5 m/s (fortissimo). A real piano spans a wider
    // range; this keeps soft notes audible, and the two-layer felt makes the tone carry the dynamics.
    const double x = std::pow (std::clamp (velocity, 1, 127) / 127.0, (double) controls.dynamicsCurve);
    return 0.8 * std::pow (5.0 / 0.8, x);
}

void PianoInstrument::noteOn (int midiNote, int velocity)
{
    if (midiNote < 0 || midiNote > 127)
        return;
    if (velocity <= 0)
    {
        noteOff (midiNote);
        return;
    }
    keyHeld[(size_t) midiNote] = true;
    piano->strike (midiNote, hammerVelocityFor (velocity));
}

void PianoInstrument::holdKeySilently (int midiNote)
{
    if (midiNote < 0 || midiNote > 127)
        return;
    keyHeld[(size_t) midiNote] = true;
    piano->setKeyDown (midiNote, true);
}

void PianoInstrument::noteOff (int midiNote)
{
    if (midiNote < 0 || midiNote > 127 || ! keyHeld[(size_t) midiNote])
        return;
    keyHeld[(size_t) midiNote] = false;
    piano->setKeyDown (midiNote, false);
}

void PianoInstrument::setSustain (bool pedalDown)
{
    piano->setSustain (pedalDown);
}

void PianoInstrument::setControls (const PianoControls& c)
{
    if (c.hammerHardness != controls.hammerHardness)
        piano->setHammerHardness (c.hammerHardness);
    if (c.damperStrength != controls.damperStrength)
        piano->setDamperStrength (c.damperStrength);
    piano->setActionNoise (c.actionNoise);
    controls = c;
}

void PianoInstrument::process (float* const* outputs, int numChannels, int numSamples)
{
    if (numChannels <= 0 || numSamples <= 0)
        return;

    for (int offset = 0; offset < numSamples;)
    {
        const int n = std::min (numSamples - offset, (int) scratch.size());
        float* left = outputs[0] + offset;
        float* right = numChannels > 1 ? outputs[1] + offset : scratch.data();
        piano->process (left, right, n);

        const float g = controls.outputGain;
        if (numChannels > 1)
            for (int i = 0; i < n; ++i) { left[i] = safetyClip (left[i] * g); right[i] = safetyClip (right[i] * g); }
        else
            for (int i = 0; i < n; ++i) left[i] = safetyClip (0.5f * g * (left[i] + right[i]));
        offset += n;
    }

    for (int ch = 2; ch < numChannels; ++ch)
        std::fill (outputs[ch], outputs[ch] + numSamples, 0.0f);

    activeStrings.store (piano->activeStringCount(), std::memory_order_relaxed);
}

} // namespace openpiano
