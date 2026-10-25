/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

// Checks the 88-key piano: every key sounds at its stretch-tuned pitch, nothing blows up
// with every key and the pedal down, dampers stop notes, and undamped strings resonate.

#include "openpiano/PianoInstrument.h"
#include "piano/ModalPiano.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace openpiano;

namespace
{
    constexpr int kRate = 48000;
    int failures = 0;

    void check (bool ok, const char* what)
    {
        if (! ok)
        {
            std::printf ("FAIL  %s\n", what);
            ++failures;
        }
    }

    // Renders n samples, mixed to mono
    std::vector<float> render (PianoInstrument& p, int n)
    {
        std::vector<float> l ((size_t) n), r ((size_t) n), out ((size_t) n);
        for (int pos = 0; pos < n; pos += 256)
        {
            const int len = std::min (256, n - pos);
            float* ch[] = { l.data() + pos, r.data() + pos };
            p.process (ch, 2, len);
        }
        for (size_t i = 0; i < out.size(); ++i)
            out[i] = 0.5f * (l[i] + r[i]);
        return out;
    }

    double rms (const std::vector<float>& x, size_t from, size_t to)
    {
        double s = 0;
        for (size_t i = from; i < to; ++i)
            s += (double) x[i] * x[i];
        return std::sqrt (s / double (std::max<size_t> (1, to - from)));
    }

    // Power of a Hann-windowed segment at one frequency (Goertzel)
    double powerAt (const std::vector<float>& x, size_t from, size_t len, double f)
    {
        const double w = 2 * M_PI * f / kRate, c = 2 * std::cos (w);
        double s1 = 0, s2 = 0;
        for (size_t i = 0; i < len; ++i)
        {
            const double win = 0.5 - 0.5 * std::cos (2 * M_PI * (double) i / (double) len);
            const double s0 = x[from + i] * win + c * s1 - s2;
            s2 = s1;
            s1 = s0;
        }
        return s1 * s1 + s2 * s2 - c * s1 * s2;
    }

    // Strongest frequency within +-30 cents of a guess, to 0.1 cent
    double measurePitch (const std::vector<float>& x, double guess)
    {
        const size_t from = (size_t) (0.3 * kRate), len = (size_t) (1.8 * kRate);
        double best = guess, bestPower = -1;
        for (double cents = -30; cents <= 30; cents += 0.1)
        {
            const double f = guess * std::pow (2.0, cents / 1200);
            const double pw = powerAt (x, from, len, f);
            if (pw > bestPower) { bestPower = pw; best = f; }
        }
        return best;
    }
}

int main()
{
    // 1. Pitch: each key against its stretch-tuned target
    double worst = 0;
    int worstNote = 0;
    {
        PianoInstrument piano;
        piano.prepare (kRate, 256);
        const auto& scale = piano.model().scale();
        for (const auto& key : scale)
        {
            piano.reset();
            piano.noteOn (key.midiNote, 90);
            const auto audio = render (piano, (int) (2.1 * kRate));
            // Bass notes: check a higher partial (near 300 Hz), which pins the pitch more finely
            const int partial = std::clamp ((int) std::round (300.0 / key.targetFrequency), 1, 8);
            const auto& s = key.strings[0];
            const double target = key.targetFrequency * partial
                                * std::sqrt ((1 + s.inharmonicity * partial * partial) / (1 + s.inharmonicity));
            const double measured = measurePitch (audio, target);
            const double cents = 1200 * std::log2 (measured / target);
            if (std::abs (cents) > std::abs (worst)) { worst = cents; worstNote = key.midiNote; }
            if (std::abs (cents) > 2.0)
                std::printf ("  MIDI %d: %+.2f cents\n", key.midiNote, cents);
        }
        std::printf ("pitch: worst error %.2f cents (MIDI %d)\n", worst, worstNote);
        check (std::abs (worst) < 2.0, "every key within 2 cents of its target pitch");
    }

    // 2. Everything at once, pedal down: stays finite and below clipping
    {
        PianoInstrument piano;
        piano.prepare (kRate, 256);
        piano.setSustain (true);
        for (int n = piano.lowestNote(); n <= piano.highestNote(); ++n)
            piano.noteOn (n, 127);
        const auto audio = render (piano, 3 * kRate);
        bool finite = true;
        float peak = 0;
        for (float v : audio) { finite = finite && std::isfinite (v); peak = std::max (peak, std::fabs (v)); }
        std::printf ("all 88 keys fortissimo with pedal: peak %.1f dBFS\n", 20 * std::log10 (peak));
        check (finite, "output stays finite with every key down");
        check (peak <= 1.0f, "output never goes past full scale");
        check (rms (audio, audio.size() - 4800, audio.size()) > 1e-4, "still ringing after 3 s with the pedal down");
    }

    // 3. Dampers: a released C4 is gone within a second; a held one is not
    {
        PianoInstrument piano;
        piano.prepare (kRate, 256);
        piano.noteOn (60, 100);
        auto a = render (piano, kRate / 2);
        piano.noteOff (60);
        auto b = render (piano, kRate);
        const double before = rms (a, a.size() - 4800, a.size()), after = rms (b, b.size() - 4800, b.size());
        std::printf ("damper: %.1f dB down 1 s after release\n", 20 * std::log10 (after / before));
        check (after < before * 1e-3, "released note decays by 60 dB within a second");
    }

    // 4. Sympathetic resonance: C3 struck short with C4 E4 G4 silently held, versus not held
    {
        double level[2];
        for (int held = 0; held < 2; ++held)
        {
            PianoInstrument piano;
            piano.prepare (kRate, 256);
            if (held)
                for (int n : { 60, 64, 67 })
                    piano.holdKeySilently (n);
            piano.noteOn (48, 110);
            render (piano, kRate / 4);
            piano.noteOff (48);
            auto tail = render (piano, kRate);
            level[held] = rms (tail, tail.size() - 9600, tail.size());
        }
        std::printf ("sympathetic: held strings ring %.1f dB above the damped case\n", 20 * std::log10 (level[1] / level[0]));
        check (level[1] > level[0] * 30, "held undamped strings pick up the struck note");
    }

    // 5. Repeatable: two renders of the same thing are identical
    {
        PianoInstrument a, b;
        a.prepare (kRate, 256);
        b.prepare (kRate, 128);
        a.noteOn (64, 80);
        b.noteOn (64, 80);
        auto x = render (a, kRate / 2);
        std::vector<float> y;
        {
            std::vector<float> l (kRate / 2), r (kRate / 2);
            for (int pos = 0; pos < kRate / 2; pos += 128)
            {
                float* ch[] = { l.data() + pos, r.data() + pos };
                b.process (ch, 2, std::min (128, kRate / 2 - pos));
            }
            for (size_t i = 0; i < l.size(); ++i)
                y.push_back (0.5f * (l[i] + r[i]));
        }
        check (x == y, "same output for any block size");
    }

    std::printf ("\n%s\n", failures == 0 ? "All piano checks passed." : "FAILED");
    return failures == 0 ? 0 : 1;
}
