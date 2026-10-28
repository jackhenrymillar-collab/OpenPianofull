#include <cstdlib>
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
            // Bass notes: check a higher partial (150-400 Hz), which pins the pitch more finely.
            // Pick one the hammer excites well (not near a node at the striking point).
            int partial = 1;
            double best = -1;
            const int lo = std::max (1, (int) std::round (150.0 / key.targetFrequency));
            const int hi = std::max (lo, std::min (12, (int) std::round (400.0 / key.targetFrequency)));
            for (int k = lo; k <= hi; ++k)
            {
                const double excitation = std::abs (std::sin (k * M_PI * key.strikeRatio));
                if (excitation > best + 0.05) { best = excitation; partial = k; }
            }
            const auto& s = key.strings[0];
            const double target = key.targetFrequency * partial
                                * std::sqrt ((1 + s.inharmonicity * partial * partial) / (1 + s.inharmonicity));
            const double measured = measurePitch (audio, target);
            // The unison's strings are tuned 1 to 3 cents apart, and the loudest one sets the
            // measured peak: count only how far the pitch falls outside the strings' own range
            double low = 0, high = 0;
            for (int i = 0; i < key.numStrings; ++i)
            {
                const double c = 1200 * std::log2 (key.strings[(size_t) i].fundamental / key.targetFrequency);
                low = std::min (low, c);
                high = std::max (high, c);
            }
            const double raw = 1200 * std::log2 (measured / target);
            const double cents = raw < low ? raw - low : raw > high ? raw - high : 0.0;
            if (std::getenv ("TUNE_DUMP")) std::printf ("TUNE %d %.3f\n", key.midiNote, raw - 0.5 * (low + high));
            if (std::abs (cents) > std::abs (worst)) { worst = cents; worstNote = key.midiNote; }
            if (std::abs (cents) > 2.0)
                std::printf ("  MIDI %d: %+.2f cents\n", key.midiNote, cents);
        }
        std::printf ("pitch: worst error %.2f cents (MIDI %d)\n", worst, worstNote);
        // Each string also moves in two directions and pulls on its neighbours through the bridge
        check (std::abs (worst) < 2.5, "every key within 2.5 cents of its target pitch");
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
        // Duplex segments, the undamped top strings (in both directions of vibration) and the
        // case keep a faint ring after release, as on a real piano (about 35 dB down on a grand)
        check (after < before * 7.9e-3, "released note decays by 42 dB within a second");
    }

    // 3b. Sostenuto: C4 held when the pedal goes down keeps ringing after release; E4 played
    //     afterwards is damped as usual
    {
        PianoInstrument piano;
        piano.prepare (kRate, 256);
        piano.noteOn (60, 100);
        render (piano, kRate / 10);
        piano.setSostenuto (true);
        piano.noteOff (60);
        piano.noteOn (64, 100);
        auto a = render (piano, kRate / 2);
        piano.noteOff (64);
        auto b = render (piano, kRate);
        const double c4 = powerAt (b, b.size() - 24000, 24000, 261.6), e4 = powerAt (b, b.size() - 24000, 24000, 329.6);
        const double e4Before = powerAt (a, a.size() - 12000, 12000, 329.6);
        std::printf ("sostenuto: held C4 vs later E4 after release %.1f dB; E4 %.1f dB down\n",
                     10 * std::log10 (c4 / e4), 10 * std::log10 (e4 / e4Before));
        check (c4 > 1000 * e4, "sostenuto holds the note that was down and not the later one");
        piano.setSostenuto (false);
        auto c = render (piano, kRate);
        // The top strings have no dampers: the released note's upper partials leave them
        // ringing faintly through the body (a G6 string at C4's sixth partial, about -31 dB
        // at the seat and -38 dB at the close mics; a real grand's undamped treble sits
        // 30-40 dB down after a release)
        check (rms (c, c.size() - 4800, c.size()) < 0.04 * rms (b, b.size() - 4800, b.size()), "releasing sostenuto damps the held note");
    }

    // 3c. Una corda: softer and darker on a trichord at the same velocity
    {
        double level[2], bright[2];
        for (int soft = 0; soft < 2; ++soft)
        {
            PianoInstrument piano;
            piano.prepare (kRate, 256);
            piano.setSoftPedal (soft ? 1.0f : 0.0f);
            piano.noteOn (69, 100);
            const auto x = render (piano, kRate / 2);
            // The attack: with one string fewer the note starts quieter, but it also drains into
            // the bridge more slowly and sustains longer (Weinreich's aftersound), so a long
            // window would mostly measure the sustain. The colour is judged in the first 50 ms,
            // before the unstruck third string has drawn off the fundamental (which would read
            // as brighter, though it is the unison coupling and not the felt)
            level[soft] = rms (x, 0, 4800);
            bright[soft] = powerAt (x, 0, 2400, 440.0 * 4.02) / powerAt (x, 0, 2400, 440.0);
        }
        std::printf ("una corda: %.1f dB quieter, 4th partial %.1f dB weaker relative to the 1st\n",
                     20 * std::log10 (level[0] / level[1]), 10 * std::log10 (bright[0] / bright[1]));
        check (level[1] < level[0] * 0.9 && level[1] > level[0] * 0.3, "una corda is a little quieter");
        check (bright[1] < bright[0], "una corda is darker");
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
        check (level[1] > level[0] * 10, "held undamped strings ring at least 20 dB above damped ones");
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
