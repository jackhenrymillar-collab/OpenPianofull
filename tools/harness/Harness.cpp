/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#include "Harness.h"

#include "openpiano/PianoInstrument.h"
#include "dr_wav.h"

#include <algorithm>
#include <cmath>

namespace openpiano::harness
{

std::vector<Scenario> standardScenarios (int sr)
{
    auto at = [sr] (double seconds) { return (int64_t) std::llround (seconds * sr); };
    std::vector<Scenario> list;

    // Single notes across the playable range, held for 3 s then released
    struct Single { const char* label; int note; int velocity; };
    const Single singles[] = {
        { "A0", 21, 80 }, { "C2", 36, 80 }, { "C3", 48, 80 }, { "C4", 60, 80 },
        { "A4", 69, 80 }, { "C5", 72, 80 }, { "C4", 60, 30 }, { "C4", 60, 127 },
    };
    for (auto& s : singles)
    {
        Scenario sc;
        sc.name = std::string ("single_") + s.label + "_v" + std::to_string (s.velocity);
        sc.description = std::string ("One ") + s.label + " at velocity " + std::to_string (s.velocity)
                       + ", held 3 s, then released";
        sc.seconds = 4.5;
        sc.events = { { at (0.01), Event::NoteOn, s.note, s.velocity }, { at (3.0), Event::NoteOff, s.note } };
        list.push_back (sc);
    }

    {
        Scenario sc { "chord_C_major", "C3 E3 G3 C4 at velocity 90, held 3 s", 4.5, {} };
        for (int n : { 48, 52, 55, 60 })
            sc.events.push_back ({ at (0.01), Event::NoteOn, n, 90 });
        for (int n : { 48, 52, 55, 60 })
            sc.events.push_back ({ at (3.0), Event::NoteOff, n });
        list.push_back (sc);
    }

    {
        Scenario sc { "repeated_C4", "C4 struck 8 times, 125 ms apart, without releasing the damper", 3.0, {} };
        for (int i = 0; i < 8; ++i)
            sc.events.push_back ({ at (0.01 + 0.125 * i), Event::NoteOn, 60, 100 });
        sc.events.push_back ({ at (2.0), Event::NoteOff, 60 });
        list.push_back (sc);
    }

    {
        Scenario sc { "pedal_arpeggio", "Sustain down, C2 G2 C3 E3 G3 C4 arpeggio with short key presses, pedal up at 4 s", 5.5, {} };
        sc.events.push_back ({ at (0.0), Event::PedalDown });
        const int notes[] = { 36, 43, 48, 52, 55, 60 };
        for (int i = 0; i < 6; ++i)
        {
            sc.events.push_back ({ at (0.01 + 0.25 * i), Event::NoteOn, notes[i], 85 });
            sc.events.push_back ({ at (0.21 + 0.25 * i), Event::NoteOff, notes[i] });
        }
        sc.events.push_back ({ at (4.0), Event::PedalUp });
        list.push_back (sc);
    }

    {
        Scenario sc { "above_range_C6", "C6, above the upstream range: currently silent", 1.0, {} };
        sc.events = { { at (0.01), Event::NoteOn, 84, 100 }, { at (0.5), Event::NoteOff, 84 } };
        list.push_back (sc);
    }

    {
        Scenario sc { "stress_all_keys", "Every playable key struck at once (CPU benchmark)", 3.0, {} };
        for (int n = 21; n <= 72; ++n)
            sc.events.push_back ({ at (0.01), Event::NoteOn, n, 100 });
        sc.benchmarkOnly = true;
        list.push_back (sc);
    }

    for (auto& sc : list)
        std::stable_sort (sc.events.begin(), sc.events.end(), [] (const Event& a, const Event& b) { return a.sample < b.sample; });

    return list;
}

static void apply (Instrument& instrument, const Event& e)
{
    switch (e.type)
    {
        case Event::NoteOn:    instrument.noteOn (e.note, e.velocity); break;
        case Event::NoteOff:   instrument.noteOff (e.note); break;
        case Event::PedalDown: instrument.setSustain (true); break;
        case Event::PedalUp:   instrument.setSustain (false); break;
    }
}

std::vector<float> renderEngine (const Scenario& sc, int sampleRate, unsigned workerThreads, int blockSize)
{
    PianoInstrument piano (workerThreads);
    piano.prepare (sampleRate, blockSize);

    const auto total = (int64_t) std::llround (sc.seconds * sampleRate);
    std::vector<float> out ((size_t) total, 0.0f);
    size_t nextEvent = 0;

    for (int64_t blockStart = 0; blockStart < total; blockStart += blockSize)
    {
        const int64_t blockEnd = std::min (total, blockStart + blockSize);
        int64_t pos = blockStart;

        // Split the block at each event, as the plugin does
        while (pos < blockEnd)
        {
            while (nextEvent < sc.events.size() && sc.events[nextEvent].sample <= pos)
                apply (piano, sc.events[nextEvent++]);

            int64_t segmentEnd = blockEnd;
            if (nextEvent < sc.events.size())
                segmentEnd = std::min (segmentEnd, sc.events[nextEvent].sample);

            float* channel = out.data() + pos;
            piano.process (&channel, 1, (int) (segmentEnd - pos));
            pos = segmentEnd;
        }
    }
    return out;
}

bool writeWav (const std::string& path, const std::vector<float>& samples, int sampleRate)
{
    drwav_data_format format {};
    format.container = drwav_container_riff;
    format.format = DR_WAVE_FORMAT_IEEE_FLOAT;
    format.channels = 1;
    format.sampleRate = (drwav_uint32) sampleRate;
    format.bitsPerSample = 32;

    drwav wav;
    if (! drwav_init_file_write (&wav, path.c_str(), &format, nullptr))
        return false;
    const auto written = drwav_write_pcm_frames (&wav, samples.size(), samples.data());
    drwav_uninit (&wav);
    return written == samples.size();
}

} // namespace openpiano::harness
