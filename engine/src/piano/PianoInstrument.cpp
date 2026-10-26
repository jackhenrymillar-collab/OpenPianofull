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
    this->sampleRate = sampleRate;
    piano::ModalPiano::Settings settings;
    settings.age = controls.age;
    settings.closeMics = closeMics;
    piano->prepare (sampleRate, settings);
    builtAge = controls.age;
    ageSettleSamples = 0;
    applyLiveControls();
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
    // Hammer speed from 0.6 m/s (pianissimo) to 5 m/s (fortissimo), which gives about the
    // loudness range of recorded grands (10 dB in the bass to 15-18 dB in the treble).
    const double x = std::pow (std::clamp (velocity, 1, 127) / 127.0, (double) controls.dynamicsCurve);
    return 0.6 * std::pow (5.0 / 0.6, x);
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
    pedal = pedalDown;
    piano->setSustain (pedalDown);
}

void PianoInstrument::setSostenuto (bool pedalDown)
{
    piano->setSostenuto (pedalDown);
}

void PianoInstrument::setSoftPedal (float amount)
{
    softAmount = amount;
    piano->setSoftPedal (amount);
}

void PianoInstrument::applyLiveControls()
{
    // An old piano's dampers are worn and let notes ring on; its action knocks and rattles more
    piano->setHammerHardness (controls.hammerHardness);
    piano->setDamperStrength (controls.damperStrength * (1.0 - 0.5 * controls.age));
    piano->setActionNoise (controls.actionNoise * (1.0 + 1.5 * controls.age));
    piano->setRoom (controls.room);
}

void PianoInstrument::setControls (const PianoControls& c)
{
    const bool changed = c.hammerHardness != controls.hammerHardness || c.damperStrength != controls.damperStrength
                         || c.actionNoise != controls.actionNoise || c.age != controls.age
                         || c.room != controls.room;
    if (c.age != controls.age)
        ageSettleSamples = (int) (0.25 * sampleRate);
    controls = c;
    if (changed)
        applyLiveControls();
}

void PianoInstrument::process (float* const* outputs, int numChannels, int numSamples)
{
    if (numChannels <= 0 || numSamples <= 0)
        return;

    if (ageSettleSamples > 0)
    {
        ageSettleSamples -= numSamples;
        if (ageSettleSamples <= 0 && controls.age != builtAge)
        {
            // Rebuild with the new age. Sounding notes stop; held keys and the pedal carry over.
            piano::ModalPiano::Settings settings;
            settings.age = controls.age;
            settings.closeMics = closeMics;
            piano->prepare (sampleRate, settings);
            builtAge = controls.age;
            applyLiveControls();
            for (int note = 0; note < 128; ++note)
                if (keyHeld[(size_t) note])
                    piano->setKeyDown (note, true);
            piano->setSustain (pedal);
            piano->setSoftPedal (softAmount);
        }
    }

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
