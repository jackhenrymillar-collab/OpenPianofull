/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include "StringModel.h"

namespace openpiano::piano
{

constexpr int kFirstMidiNote = 21; // A0
constexpr int kNumKeys = 88;       // A0..C8

// The upstream parameter set: the same placeholder values for every string, with only the
// frequency, length and (from C5 up) the cross-section changing. Upstream only enables
// A0..C5, because the finite-difference grid gets too coarse above that at 48 kHz.
// Finding real per-string values is part of the next step of the plan.
constexpr int kLastEnabledMidiNote = 72; // C5

struct KeyParams
{
    StringParams string;
    HammerParams hammer;
};

inline KeyParams upstreamKeyParams (int midiNote)
{
    static constexpr double f0[kNumKeys] = {
        27.5, 29.14, 30.87,
        32.7, 34.65, 36.71, 38.89, 41.20, 43.65, 46.25, 49.00, 51.91, 55.00, 58.27, 61.74,
        65.41, 69.30, 73.42, 77.78, 82.41, 87.31, 92.50, 98.00, 103.83, 110.00, 116.54, 123.47,
        130.81, 138.59, 146.83, 155.56, 164.81, 174.61, 185.00, 196, 207.65, 220, 233.08, 246.94,
        261.63, 277.18, 293.66, 311.13, 329.63, 349.23, 369.99, 392.00, 415.30, 440.00, 466.16, 493.88,
        523.25, 554.37, 587.33, 622.25, 659.26, 698.46, 739.99, 783.99, 830.61, 880.00, 932.33, 987.77,
        1046.50, 1108.73, 1174.66, 1244.51, 1318.51, 1396.91, 1479.98, 1567.98, 1661.22, 1760.00, 1864.66, 1975.53,
        2093.00, 2217.46, 2349.32, 2489.02, 2637.02, 2793.83, 2959.96, 3135.96, 3322.44, 3520.00, 3729.31, 3951.07,
        4186.01
    };

    const int key = midiNote - kFirstMidiNote;
    const bool bass = midiNote < 48;         // below C3
    const double L = bass ? 1.92 : 0.96;
    const double S = midiNote < 72 ? 0.001   // below C5
                   : midiNote < 84 ? 0.0008  // C5..B5
                   : 0.0005;                 // C6 and up

    KeyParams p;
    p.string = { f0[key], L, 0.0182, S, 9e7, 0.003, 6.25e-9 };
    p.hammer = { 4.9e-03, 2.3, 1e-04, 4e08, 0.12, 0.05 };
    return p;
}

} // namespace openpiano::piano
