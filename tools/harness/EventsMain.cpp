/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

// openpiano-events: renders a text list of piano events to a stereo WAV file.
//
//   openpiano-events EVENTS OUT.wav SECONDS [AGE]
//
// Each line of EVENTS is "<seconds> on <note> <velocity>", "<seconds> off <note> 0" or
// "<seconds> ped <0|1> 0", "<seconds> sost <0|1> 0" or "<seconds> soft <0-127> 0" (sostenuto
// and una corda), in time order. Rendered at 48 kHz in 256-sample blocks with the
// default controls (Age from the optional last argument).

#include "Harness.h"
#include "openpiano/PianoInstrument.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

int main (int argc, char** argv)
{
    if (argc < 4)
    {
        std::puts ("usage: openpiano-events EVENTS OUT.wav SECONDS [AGE]");
        return 2;
    }
    constexpr int rate = 48000, block = 256;
    const double seconds = std::atof (argv[3]);

    struct Event { double t; std::string kind; int a, b; };
    std::vector<Event> events;
    FILE* f = std::fopen (argv[1], "r");
    if (f == nullptr)
    {
        std::fprintf (stderr, "cannot open %s\n", argv[1]);
        return 1;
    }
    char kind[8];
    double t;
    int a, b;
    while (std::fscanf (f, "%lf %7s %d %d", &t, kind, &a, &b) == 4)
        events.push_back ({ t, kind, a, b });
    std::fclose (f);

    openpiano::PianoInstrument piano;
    openpiano::PianoControls controls;
    controls.age = argc > 4 ? (float) std::atof (argv[4]) : 0.0f;
    piano.setControls (controls);
    piano.prepare (rate, block);

    const int total = (int) (seconds * rate);
    std::vector<float> left ((size_t) total), right ((size_t) total), out (2 * (size_t) total);
    size_t next = 0;
    for (int pos = 0; pos < total; pos += block)
    {
        while (next < events.size() && events[next].t * rate < pos + block)
        {
            const auto& e = events[next++];
            if (e.kind == "on") piano.noteOn (e.a, e.b);
            else if (e.kind == "off") piano.noteOff (e.a);
            else if (e.kind == "ped") piano.setSustain (e.a != 0);
            else if (e.kind == "sost") piano.setSostenuto (e.a != 0);
            else if (e.kind == "soft") piano.setSoftPedal ((float) e.a / 127.0f);
        }
        const int n = std::min (block, total - pos);
        float* ch[] = { left.data() + pos, right.data() + pos };
        piano.process (ch, 2, n);
    }
    for (size_t i = 0; i < (size_t) total; ++i)
    {
        out[2 * i] = left[i];
        out[2 * i + 1] = right[i];
    }
    return openpiano::harness::writeWav (argv[2], out, rate, 2) ? 0 : 1;
}
