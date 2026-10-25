/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#include "Harness.h"

#include "openpiano/LegacyPiano.h"
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
        { "A4", 69, 80 }, { "C5", 72, 80 }, { "C4", 60, 30 }, { "C4", 60, 127 }, { "C2", 36, 127 },
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

    // The rest use the full 88-key range of the new piano
    for (auto& s : { Single { "C6", 84, 80 }, Single { "C7", 96, 80 }, Single { "C8", 108, 80 } })
    {
        Scenario sc;
        sc.name = std::string ("single_") + s.label + "_v" + std::to_string (s.velocity);
        sc.description = std::string ("One ") + s.label + " at velocity " + std::to_string (s.velocity)
                       + ", held 3 s, then released";
        sc.seconds = 4.5;
        sc.events = { { at (0.01), Event::NoteOn, s.note, s.velocity }, { at (3.0), Event::NoteOff, s.note } };
        sc.fullRangeOnly = true;
        list.push_back (sc);
    }

    {
        Scenario sc { "sympathetic_pedal", "C3 struck short, once with the sustain pedal up and once with it down", 7.0, {} };
        sc.events.push_back ({ at (0.01), Event::NoteOn, 48, 100 });
        sc.events.push_back ({ at (0.30), Event::NoteOff, 48 });
        sc.events.push_back ({ at (3.0), Event::PedalDown });
        sc.events.push_back ({ at (3.51), Event::NoteOn, 48, 100 });
        sc.events.push_back ({ at (3.80), Event::NoteOff, 48 });
        sc.events.push_back ({ at (6.5), Event::PedalUp });
        sc.fullRangeOnly = true;
        list.push_back (sc);
    }

    {
        Scenario sc { "sympathetic_held_chord", "C4 E4 G4 held down silently, then C3 struck short: the held strings ring on", 4.0, {} };
        for (int n : { 60, 64, 67 })
            sc.events.push_back ({ at (0.0), Event::SilentKeyDown, n });
        sc.events.push_back ({ at (0.5), Event::NoteOn, 48, 110 });
        sc.events.push_back ({ at (0.75), Event::NoteOff, 48 });
        sc.fullRangeOnly = true;
        list.push_back (sc);
    }

    {
        Scenario sc { "scale_full_range", "A0 to C8 in octaves, then a high flourish without dampers", 9.0, {} };
        int i = 0;
        for (int n = 21; n <= 108; n += 12, ++i)
        {
            sc.events.push_back ({ at (0.01 + 0.6 * i), Event::NoteOn, n, 90 });
            sc.events.push_back ({ at (0.5 + 0.6 * i), Event::NoteOff, n });
        }
        for (int j = 0; j < 8; ++j)
        {
            const int n = 96 + (j * 5) % 12;
            sc.events.push_back ({ at (5.5 + 0.1 * j), Event::NoteOn, n, 70 });
            sc.events.push_back ({ at (5.58 + 0.1 * j), Event::NoteOff, n });
        }
        sc.fullRangeOnly = true;
        list.push_back (sc);
    }

    // The pedalled arpeggio at four touches, for comparing dynamics
    for (int velocity : { 30, 60, 90, 127 })
    {
        Scenario sc { "pedal_arpeggio_v" + std::to_string (velocity),
                      "Sustain down, C2 G2 C3 E3 G3 C4 arpeggio at velocity " + std::to_string (velocity) + ", pedal up at 4 s",
                      5.5, {} };
        sc.events.push_back ({ at (0.0), Event::PedalDown });
        const int notes[] = { 36, 43, 48, 52, 55, 60 };
        for (int i = 0; i < 6; ++i)
        {
            sc.events.push_back ({ at (0.01 + 0.25 * i), Event::NoteOn, notes[i], velocity });
            sc.events.push_back ({ at (0.21 + 0.25 * i), Event::NoteOff, notes[i] });
        }
        sc.events.push_back ({ at (4.0), Event::PedalUp });
        sc.fullRangeOnly = true;
        list.push_back (sc);
    }

    {
        Scenario sc { "stress_all_88", "All 88 keys struck at once with the pedal down (CPU benchmark)", 4.0, {} };
        sc.events.push_back ({ at (0.0), Event::PedalDown });
        for (int n = 21; n <= 108; ++n)
            sc.events.push_back ({ at (0.01), Event::NoteOn, n, 100 });
        sc.benchmarkOnly = true;
        sc.fullRangeOnly = true;
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
        case Event::SilentKeyDown:
            // Pressing a key slowly enough that the hammer never reaches the strings
            if (auto* p = dynamic_cast<PianoInstrument*> (&instrument))
                p->holdKeySilently (e.note);
            break;
    }
}

template <typename Engine, typename Process>
static void drive (Engine& piano, const Scenario& sc, int sampleRate, int blockSize, Process&& process)
{
    const auto total = (int64_t) std::llround (sc.seconds * sampleRate);
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

            process (pos, (int) (segmentEnd - pos));
            pos = segmentEnd;
        }
    }
}

std::vector<float> renderLegacy (const Scenario& sc, int sampleRate, unsigned workerThreads, int blockSize)
{
    LegacyPiano piano (workerThreads);
    piano.prepare (sampleRate, blockSize);

    std::vector<float> out ((size_t) std::llround (sc.seconds * sampleRate), 0.0f);
    drive (piano, sc, sampleRate, blockSize, [&] (int64_t pos, int n)
    {
        float* channel = out.data() + pos;
        piano.process (&channel, 1, n);
    });
    return out;
}

std::vector<float> renderPiano (const Scenario& sc, int sampleRate, int blockSize)
{
    PianoInstrument piano;
    piano.prepare (sampleRate, blockSize);

    const auto total = (size_t) std::llround (sc.seconds * sampleRate);
    std::vector<float> left (total), right (total), out (2 * total);
    drive (piano, sc, sampleRate, blockSize, [&] (int64_t pos, int n)
    {
        float* channels[] = { left.data() + pos, right.data() + pos };
        piano.process (channels, 2, n);
    });
    for (size_t i = 0; i < total; ++i)
    {
        out[2 * i] = left[i];
        out[2 * i + 1] = right[i];
    }
    return out;
}

bool writeWav (const std::string& path, const std::vector<float>& samples, int sampleRate, int channels)
{
    drwav_data_format format {};
    format.container = drwav_container_riff;
    format.format = DR_WAVE_FORMAT_IEEE_FLOAT;
    format.channels = (drwav_uint32) channels;
    format.sampleRate = (drwav_uint32) sampleRate;
    format.bitsPerSample = 32;

    drwav wav;
    if (! drwav_init_file_write (&wav, path.c_str(), &format, nullptr))
        return false;
    const auto written = drwav_write_pcm_frames (&wav, samples.size() / (size_t) channels, samples.data());
    drwav_uninit (&wav);
    return written == samples.size() / (size_t) channels;
}

} // namespace openpiano::harness
