/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

namespace openpiano
{

// Common interface for every instrument mode (Piano now, Strings and Rhodes later).
//
// Threading contract:
//  - prepare() may allocate and is called from a non-realtime thread.
//  - All other methods are called from the audio thread only and must not allocate or lock.
//  - Events (noteOn, noteOff, the pedals) take effect from the next sample rendered by process().
//    Callers that want sample-accurate timing split their block at each event.
class Instrument
{
public:
    virtual ~Instrument() = default;

    virtual const char* name() const = 0;

    virtual void prepare (double sampleRate, int maxBlockSize) = 0;

    // Silences everything immediately (e.g. on transport reset).
    virtual void reset() = 0;

    // velocity is the MIDI velocity, 1..127.
    virtual void noteOn (int midiNote, int velocity) = 0;
    virtual void noteOff (int midiNote) = 0;
    // releaseVelocity is the MIDI note-off velocity, 1..127 (64 when the keyboard sends none):
    // how fast the key was let go. Instruments that don't use it fall back to a plain note-off.
    virtual void noteOff (int midiNote, int /*releaseVelocity*/) { noteOff (midiNote); }
    virtual void setSustain (bool pedalDown) = 0;
    // Sustain pedal position, 0 = up, 1 = fully down (half pedalling). Instruments without a
    // continuous pedal treat anything past halfway as down.
    virtual void setSustainLevel (float level) { setSustain (level >= 0.5f); }
    // Middle and left pedals. Instruments without them ignore these.
    virtual void setSostenuto (bool /*pedalDown*/) {}
    virtual void setSoftPedal (float /*amount*/) {} // 0 = up, 1 = fully down

    // Overwrites numSamples samples in each output channel. numSamples may be any length,
    // including more than the maxBlockSize given to prepare().
    virtual void process (float* const* outputs, int numChannels, int numSamples) = 0;

    // Lowest and highest MIDI notes this instrument plays.
    virtual int lowestNote() const = 0;
    virtual int highestNote() const = 0;

    // Number of vibrating elements currently being simulated (for the UI).
    virtual int activeVoiceCount() const = 0;
};

} // namespace openpiano
