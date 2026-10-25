/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#include "openpiano/PianoInstrument.h"

#include "NoteTable.h"
#include "StringModel.h"

#include <algorithm>
#include <cmath>
#include <thread>

namespace openpiano
{

using namespace piano;

unsigned PianoInstrument::defaultWorkerThreads()
{
    // Leave a core for the host. More than 7 workers doesn't help with this many strings.
    const unsigned hw = std::thread::hardware_concurrency();
    return hw <= 2 ? 0u : std::min (hw - 2, 7u);
}

PianoInstrument::PianoInstrument (unsigned numWorkerThreads)
    : pool (std::make_unique<ParallelFor> (numWorkerThreads))
{
}

PianoInstrument::~PianoInstrument() = default;

int PianoInstrument::lowestNote() const { return kFirstMidiNote; }
int PianoInstrument::highestNote() const { return kLastEnabledMidiNote; }

void PianoInstrument::prepare (double newSampleRate, int maxBlockSize)
{
    sampleRate = newSampleRate;
    chunkSize = std::max (32, maxBlockSize);

    strings.clear();
    for (int note = lowestNote(); note <= highestNote(); ++note)
    {
        const auto p = upstreamKeyParams (note);
        strings.push_back (std::make_unique<StringModel> ((int) sampleRate, p.string, p.hammer));
        strings.back()->setHammerHardness (controls.hammerHardness);
    }

    stringBuffers.assign (strings.size(), std::vector<double> ((size_t) chunkSize, 0.0));
    keyHeld.fill (false);
    sustainDown = false;
    activeVoices = 0;
}

void PianoInstrument::reset()
{
    for (auto& s : strings)
        s->reset();
    keyHeld.fill (false);
    sustainDown = false;
}

StringModel* PianoInstrument::stringFor (int midiNote)
{
    if (midiNote < lowestNote() || midiNote > highestNote() || strings.empty())
        return nullptr;
    return strings[(size_t) (midiNote - lowestNote())].get();
}

double PianoInstrument::hammerVelocityFor (int velocity) const
{
    velocity = std::clamp (velocity, 1, 127);
    if (controls.dynamicsCurve == 1.0f)
        return velocity / 30.0; // upstream mapping: velocity 127 -> 4.2 m/s
    return (127.0 / 30.0) * std::pow (velocity / 127.0, (double) controls.dynamicsCurve);
}

void PianoInstrument::noteOn (int midiNote, int velocity)
{
    if (midiNote < 0 || midiNote > 127)
        return;
    keyHeld[(size_t) midiNote] = true;
    if (auto* s = stringFor (midiNote))
        s->hit (hammerVelocityFor (velocity));
}

void PianoInstrument::noteOff (int midiNote)
{
    if (midiNote < 0 || midiNote > 127)
        return;
    keyHeld[(size_t) midiNote] = false;
    if (sustainDown)
        return;
    if (auto* s = stringFor (midiNote))
        s->damp (controls.damperStrength);
}

void PianoInstrument::setSustain (bool pedalDown)
{
    const bool released = sustainDown && ! pedalDown;
    sustainDown = pedalDown;
    if (! released)
        return;

    // Pedal up: damp every string whose key isn't held
    for (int note = lowestNote(); note <= highestNote(); ++note)
        if (! keyHeld[(size_t) note])
            stringFor (note)->damp (controls.damperStrength);
}

void PianoInstrument::setControls (const PianoControls& c)
{
    if (c.hammerHardness != controls.hammerHardness)
        for (auto& s : strings)
            s->setHammerHardness (c.hammerHardness);
    controls = c;
}

void PianoInstrument::process (float* const* outputs, int numChannels, int numSamples)
{
    if (numChannels <= 0)
        return;

    for (int offset = 0; offset < numSamples; offset += chunkSize)
        renderChunk (outputs[0] + offset, std::min (chunkSize, numSamples - offset));

    for (int ch = 1; ch < numChannels; ++ch)
        std::copy (outputs[0], outputs[0] + numSamples, outputs[ch]);

    int voices = 0;
    for (auto& s : strings)
        voices += s->isActive() ? 1 : 0;
    activeVoices.store (voices, std::memory_order_relaxed);
}

void PianoInstrument::renderChunk (float* out, int numSamples)
{
    if (strings.empty())
    {
        std::fill (out, out + numSamples, 0.0f);
        return;
    }

    auto job = [this, numSamples] (int stringIndex)
    {
        auto& s = *strings[(size_t) stringIndex];
        double* buffer = stringBuffers[(size_t) stringIndex].data();
        for (int i = 0; i < numSamples; ++i)
            buffer[i] = s.nextSample();
    };
    pool->run ((int) strings.size(), job);

    // Mix in string order, accumulating in float exactly as upstream did
    const float gain = kEngineGain * controls.outputGain;
    for (int i = 0; i < numSamples; ++i)
    {
        float sample = 0;
        for (auto& buffer : stringBuffers)
            sample += buffer[(size_t) i];
        out[i] = gain * sample;
    }
}

} // namespace openpiano
