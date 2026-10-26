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
//  - Events (noteOn, noteOff, setSustain) take effect from the next sample rendered by process().
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
    virtual void setSustain (bool pedalDown) = 0;

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
