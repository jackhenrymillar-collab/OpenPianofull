/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

// Drives the original OpenPiano code (frozen copy in tools/upstream) through a scenario.
// The event handling mirrors upstream's PluginProcessor: velocity/30 m/s on note-on,
// damp on note-off unless the pedal is down, damp unheld strings when the pedal comes up.

#include "Harness.h"
#include "piano.h" // upstream

#include <array>
#include <cmath>

namespace openpiano::harness
{

std::vector<float> renderUpstream (const Scenario& sc, int sampleRate)
{
    Piano piano (sampleRate, 256, 1);

    const auto total = (int64_t) std::llround (sc.seconds * sampleRate);
    std::vector<float> out ((size_t) total, 0.0f);
    std::array<bool, 128> held {};
    bool pedal = false;
    size_t nextEvent = 0;

    auto stringFor = [&] (int note) -> PianoString*
    {
        const int idx = note - MIDI_NOTE_OFFSET;
        return (idx >= 0 && idx < N_STRINGS) ? piano.strings[idx] : nullptr;
    };

    for (int64_t n = 0; n < total; ++n)
    {
        while (nextEvent < sc.events.size() && sc.events[nextEvent].sample <= n)
        {
            const auto& e = sc.events[nextEvent++];
            switch (e.type)
            {
                case Event::NoteOn:
                    held[(size_t) e.note] = true;
                    if (auto* s = stringFor (e.note))
                        s->hit (e.velocity / 30.0);
                    break;
                case Event::NoteOff:
                    held[(size_t) e.note] = false;
                    if (! pedal)
                        if (auto* s = stringFor (e.note))
                            s->damp();
                    break;
                case Event::PedalDown:
                    pedal = true;
                    break;
                case Event::PedalUp:
                    if (pedal)
                        for (int i = 0; i < N_STRINGS; ++i)
                            if (! held[(size_t) (i + MIDI_NOTE_OFFSET)])
                                piano.strings[i]->damp();
                    pedal = false;
                    break;
            }
        }
        out[(size_t) n] = 150.0f * piano.get_next_sample (1.0f);
    }
    return out;
}

} // namespace openpiano::harness
