/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include <array>
#include <vector>

namespace openpiano::piano
{

// The physical design of one key: its strings, hammer and tuning.
struct StringSpec
{
    double length;        // Speaking length [m]
    double coreDiameter;  // Steel core (or the whole wire, if plain) [m]
    double outerDiameter; // Including copper winding [m]
    double linearDensity; // Mass per metre [kg/m]
    double tension;       // [N]
    double inharmonicity; // B: partial k sits at k * f1 * sqrt((1 + B k^2) / (1 + B))
    double fundamental;   // Frequency of the first partial [Hz]
    bool wound;
};

struct KeySpec
{
    int midiNote;
    int numStrings;                  // 1, 2 or 3
    std::array<StringSpec, 3> strings;
    double targetFrequency;          // Tuned first partial, before unison detuning [Hz]
    double stretchCents;             // Tuning relative to equal temperament (A4 = 440 Hz)

    double hammerMass;               // [kg]
    double outerFeltDepth;           // [m] soft outer felt before the hard core
    double feltStiffness;            // Q0 in F = Q0 * (compression in mm)^p  [N]
    double feltExponent;             // p
    double feltDamping;              // [s] the felt's loss: F = Q0 (x^p + feltDamping * d(x^p)/dt)
    double strikeRatio;              // Striking point as a fraction of the speaking length
    double hammerWidth;              // Contact width along the string [m]

    bool hasDamper;                  // The top notes of a piano have no dampers
    bool bassBridge;                 // Which bridge the strings cross
};

// Builds the whole 88-key scale. See ScaleDesign.cpp for sources.
//   unisonSpreadCents: how far apart the strings of one key are tuned (0 = identical)
std::vector<KeySpec> buildConcertGrandScale (double unisonSpreadCents);

constexpr int kLowestKey = 21;  // A0
constexpr int kHighestKey = 108; // C8

} // namespace openpiano::piano
