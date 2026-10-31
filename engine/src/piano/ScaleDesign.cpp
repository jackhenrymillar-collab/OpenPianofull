/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#include "ScaleDesign.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdint>

namespace openpiano::piano
{

// Sources
// -------
// Scale: the Steinway D of Chabassier & Duruflé (INRIA RT-0425, 2012), appendix A: speaking lengths,
//   wire gauges, tensions and striking points of every key. The jump at key 21 is the break between
//   the bass bridge and the treble bridge.
// String layout: monochords A0-E1, bichords F1-E2 (bass bridge, copper-wound), plain trichords F2-C8. Dampers stop at F6, as on concert grands.
// Hammers: the same report's per-key head mass, felt exponent and stiffness laws (see below).
// Wire: steel E = 2.0e11 Pa, density 7850 kg/m^3; copper winding 8960 kg/m^3. Inharmonicity
//   B = pi^3 E d^4 / (64 T L^2) from the steel core, the standard result for stiff strings (Fletcher 1964).
// Tuning: A4 = 440 Hz. A3-A4 is the temperament octave; every other octave is tuned
//   upwards and downwards as a beatless octave on the actual (inharmonic) partials, as an
//   aural tuner would: 4:2 in the middle, 6:3 and 8:4 in the bass, 2:1 in the top octaves.
//   This produces the stretched "Railsback" curve of a real piano.

namespace
{
    // Steinway D scale, every key from Chabassier & Durufle (Inria RT-425, 2012), appendix A:
    // speaking length, wire diameter (the steel core of the wound strings), tension at pitch and
    // the hammer's striking distance from the agraffe. A0-B0 and C8 are extended from their neighbours.
    constexpr double kLength[88] = {
        2.010, 2.009, 2.008, 2.007, 1.997, 1.981, 1.965, 1.938, 1.911, 1.879, 1.842, 1.805,
        1.762, 1.709, 1.655, 1.602, 1.548, 1.495, 1.442, 1.378, 1.837, 1.757, 1.660, 1.591,
        1.482, 1.403, 1.329, 1.259, 1.192, 1.129, 1.070, 1.013, 0.960, 0.909, 0.861, 0.816,
        0.773, 0.732, 0.694, 0.657, 0.622, 0.590, 0.559, 0.529, 0.501, 0.475, 0.450, 0.426,
        0.404, 0.383, 0.363, 0.344, 0.326, 0.308, 0.292, 0.277, 0.262, 0.249, 0.236, 0.223,
        0.211, 0.200, 0.190, 0.180, 0.171, 0.162, 0.153, 0.145, 0.138, 0.130, 0.124, 0.117,
        0.111, 0.105, 0.100, 0.095, 0.090, 0.085, 0.081, 0.076, 0.072, 0.069, 0.065, 0.062,
        0.058, 0.055, 0.052, 0.049 };
    constexpr double kDiameterMm[88] = {
        1.480, 1.480, 1.480, 1.480, 1.505, 1.506, 1.492, 1.460, 1.419, 1.370, 1.316, 1.262,
        1.207, 1.148, 1.096, 1.051, 1.012, 0.982, 0.960, 0.937, 1.117, 1.108, 1.091, 1.095,
        1.071, 1.067, 1.065, 1.063, 1.061, 1.059, 1.057, 1.053, 1.049, 1.045, 1.039, 1.033,
        1.027, 1.020, 1.013, 1.006, 0.999, 0.991, 0.984, 0.977, 0.970, 0.964, 0.958, 0.952,
        0.947, 0.941, 0.937, 0.932, 0.928, 0.924, 0.921, 0.917, 0.914, 0.910, 0.907, 0.904,
        0.901, 0.898, 0.894, 0.891, 0.887, 0.883, 0.879, 0.875, 0.870, 0.866, 0.860, 0.855,
        0.850, 0.844, 0.837, 0.831, 0.824, 0.817, 0.810, 0.802, 0.795, 0.787, 0.778, 0.770,
        0.761, 0.752, 0.743, 0.770 };
    constexpr double kTension[88] = {
        1722, 1722, 1722, 1722, 1788, 1798, 1773, 1703, 1618, 1515, 1404, 1297,
        1190, 1083, 991, 915, 853, 807, 774, 741, 799, 791, 773, 783,
        754, 754, 756, 759, 762, 764, 766, 767, 766, 765, 763, 759,
        755, 751, 746, 741, 735, 730, 725, 720, 715, 711, 707, 704,
        701, 699, 697, 696, 695, 694, 694, 694, 694, 695, 695, 696,
        697, 697, 697, 697, 697, 697, 696, 695, 693, 691, 689, 686,
        683, 679, 674, 670, 664, 659, 652, 646, 639, 631, 623, 615,
        606, 598, 588, 580 };
    // The technician's voicing of the top two octaves (C6 to C8, one entry per key). Up there
    // the hammer's blow lasts three to six periods of the string, so the fundamental sits near
    // a null of the force pulse's spectrum and a tenth more or less felt stiffness moves a
    // note's loudness by 10 dB; every concert grand is voiced hammer by hammer for that reason.
    // These factors on the felt stiffness were set once for even mezzo-forte loudness against
    // each key's neighbours (a five-key running median of the unvoiced levels, so the register
    // trend is left to the physics), within what needling and hardening can do (1/3 to 3 times).
    // Unevenness 3.6 -> 1.2 dB rms; D#6 and G#7 stay 5 and 3 dB quiet, which is not the hammer.
    constexpr double kTopVoicing[25] = {
        0.33, 0.86, 1.00, 3.00, 0.78, 0.40, 1.06, 1.00, 1.51, 0.71, 3.00, 3.00, 1.00,
        0.33, 0.53, 3.00, 0.87, 3.00, 0.55, 0.89, 1.59, 1.00, 1.55, 0.56, 1.00 };

    constexpr double kStrikeDistance[88] = {
        0.241, 0.241, 0.241, 0.241, 0.240, 0.238, 0.236, 0.233, 0.229, 0.225, 0.221, 0.217,
        0.211, 0.205, 0.199, 0.192, 0.186, 0.179, 0.173, 0.165, 0.220, 0.211, 0.199, 0.191,
        0.178, 0.168, 0.160, 0.151, 0.143, 0.136, 0.128, 0.122, 0.115, 0.109, 0.103, 0.098,
        0.093, 0.088, 0.083, 0.079, 0.075, 0.071, 0.067, 0.064, 0.060, 0.057, 0.054, 0.051,
        0.048, 0.046, 0.044, 0.041, 0.039, 0.037, 0.035, 0.033, 0.031, 0.030, 0.028, 0.027,
        0.025, 0.024, 0.023, 0.022, 0.020, 0.019, 0.018, 0.017, 0.017, 0.016, 0.015, 0.014,
        0.013, 0.013, 0.012, 0.011, 0.011, 0.010, 0.010, 0.009, 0.009, 0.008, 0.008, 0.007,
        0.007, 0.007, 0.006, 0.006 };

    double lengthOf (int key) { return kLength[key - 1]; }

    constexpr double kSteelE = 2.0e11;
    constexpr double kSteelDensity = 7850.0;
    constexpr double kCopperDensity = 8960.0;
    constexpr double kWindingPacking = 0.85;

    constexpr int kLastMonochordKey = 8;   // E1
    constexpr int kLastBassBridgeKey = 20; // E2
    constexpr int kLastWoundKey = 20;      // E2: on the D the treble bridge starts with plain wire
    constexpr int kLastDamperMidi = 89;    // F6

    double lerp (double a, double b, double t) { return a + (b - a) * t; }

    double equalTempered (int midi) { return 440.0 * std::pow (2.0, (midi - 69) / 12.0); }

    double inharmonicityOf (double coreDiameter, double tension, double length)
    {
        return std::pow (M_PI, 3) * kSteelE * std::pow (coreDiameter, 4) / (64.0 * tension * length * length);
    }

    // Frequency of partial n relative to the first partial, on a stiff string
    double partial (int n, double B) { return n * std::sqrt ((1.0 + B * n * n) / (1.0 + B)); }

    // Octave tuned so that partial 2m of the lower string beats with partial m of the upper
    // one: f_upper / f_lower. Tuners listen to higher coincident partials lower down.
    double octaveRatio (int m, double bLower, double bUpper) { return partial (2 * m, bLower) / partial (m, bUpper); }

    // Which coincident partials a tuner checks for the octave above this key
    int octaveCheck (int key)
    {
        if (key >= 37) return 2;  // 4:2 from A3 up
        if (key >= 16) return 3;  // 6:3 in the tenor and upper bass (C2 to G#3)
        return 4;                 // 8:4 in the low bass
    }

    // A deterministic per-key variation in [-1, 1], so every build sounds the same
    double keyNoise (int key, int salt)
    {
        uint32_t h = (uint32_t) key * 2654435761u ^ (uint32_t) salt * 40503u;
        h ^= h >> 15; h *= 2246822519u; h ^= h >> 13;
        return (h & 0xffff) / 32767.5 - 1.0;
    }

    struct Design
    {
        bool wound;
        double coreDiameter;
        double plainTargetTension;
        double woundTension;
    };

    Design designFor (int key, double previousPlainDiameter, double frequency, double length)
    {
        (void) previousPlainDiameter; (void) frequency; (void) length;
        Design d {};
        d.wound = key <= kLastWoundKey;
        d.coreDiameter = kDiameterMm[key - 1] * 1e-3;
        d.woundTension = kTension[key - 1];
        d.plainTargetTension = kTension[key - 1];
        return d;
    }

    // Fills in a string for a given first-partial frequency
    StringSpec makeString (const Design& d, double length, double fundamental)
    {
        StringSpec s {};
        s.length = length;
        s.wound = d.wound;
        s.coreDiameter = d.coreDiameter;
        s.fundamental = fundamental;

        if (d.wound)
        {
            s.tension = d.woundTension;
            // B depends only on core, tension and length; the flexible-string frequency is f1 / sqrt(1 + B)
            s.inharmonicity = inharmonicityOf (s.coreDiameter, s.tension, length);
            const double fFlex = fundamental / std::sqrt (1.0 + s.inharmonicity);
            s.linearDensity = s.tension / std::pow (2.0 * length * fFlex, 2);
            const double coreMu = kSteelDensity * M_PI * s.coreDiameter * s.coreDiameter / 4.0;
            const double windingArea = std::max (0.0, s.linearDensity - coreMu) / (kCopperDensity * kWindingPacking);
            s.outerDiameter = std::sqrt (s.coreDiameter * s.coreDiameter + 4.0 * windingArea / M_PI);
        }
        else
        {
            s.linearDensity = kSteelDensity * M_PI * s.coreDiameter * s.coreDiameter / 4.0;
            s.outerDiameter = s.coreDiameter;
            // Tension and B depend on each other: iterate (converges in a few steps)
            s.inharmonicity = 0.0;
            for (int i = 0; i < 6; ++i)
            {
                const double fFlex = fundamental / std::sqrt (1.0 + s.inharmonicity);
                s.tension = s.linearDensity * std::pow (2.0 * length * fFlex, 2);
                s.inharmonicity = inharmonicityOf (s.coreDiameter, s.tension, length);
            }
        }
        return s;
    }
}

std::vector<KeySpec> buildConcertGrandScale (double unisonSpreadCents)
{
    std::vector<KeySpec> keys (88);
    std::vector<Design> designs (88);

    // 1. Wire choices, made at equal temperament
    double previousPlain = 0;
    for (int k = 1; k <= 88; ++k)
    {
        const int midi = k + 20;
        designs[(size_t) k - 1] = designFor (k, previousPlain, equalTempered (midi), lengthOf (k));
        if (! designs[(size_t) k - 1].wound)
            previousPlain = designs[(size_t) k - 1].coreDiameter;
    }

    // 2. Stretch tuning on the real partials. B barely depends on the small tuning offsets,
    //    so two passes are plenty.
    std::vector<double> F (88), B (88);
    for (int k = 1; k <= 88; ++k)
    {
        F[(size_t) k - 1] = equalTempered (k + 20);
        B[(size_t) k - 1] = makeString (designs[(size_t) k - 1], lengthOf (k), F[(size_t) k - 1]).inharmonicity;
    }

    const int a3 = 57 - 20, a4 = 69 - 20; // key numbers
    for (int pass = 0; pass < 2; ++pass)
    {
        auto f = [&] (int key) -> double& { return F[(size_t) key - 1]; };
        auto b = [&] (int key) { return B[(size_t) key - 1]; };

        // Temperament octave A3..A4: a 4:2 octave divided into equal semitones
        f (a4) = 440.0;
        f (a3) = 440.0 / octaveRatio (2, b (a3), b (a4));
        for (int key = a3 + 1; key < a4; ++key)
            f (key) = f (a3) * std::pow (f (a4) / f (a3), (key - a3) / 12.0);

        // Upwards: 4:2 octaves all the way, which stretches the top like the recorded grand
        // (about +22 cents at C7, +38 at A7)
        for (int key = a4 + 1; key <= 88; ++key)
            f (key) = f (key - 12) * octaveRatio (2, b (key - 12), b (key));
        // Downwards: 4:2, then 6:3 and 8:4 octaves as the strings get lower
        for (int key = a3 - 1; key >= 1; --key)
            f (key) = f (key + 12) / octaveRatio (octaveCheck (key), b (key), b (key + 12));

        for (int key = 1; key <= 88; ++key)
            B[(size_t) key - 1] = makeString (designs[(size_t) key - 1], lengthOf (key), f (key)).inharmonicity;
    }

    // 3. Keys: strings (with unison spread), hammers, dampers
    for (int k = 1; k <= 88; ++k)
    {
        auto& key = keys[(size_t) k - 1];
        key.midiNote = k + 20;
        key.numStrings = k <= kLastMonochordKey ? 1 : k <= kLastBassBridgeKey ? 2 : 3;
        key.targetFrequency = F[(size_t) k - 1];
        key.stretchCents = 1200.0 * std::log2 (key.targetFrequency / equalTempered (key.midiNote));

        // Unison: trichord offsets around the centre; the exact amount varies per key.
        // One tuner's rule for the whole keyboard: a careful concert tuning leaves the strings of a
        // unison a few tenths of a cent either side of the centre. That close, the bridge locks
        // them together: the in-step motion drains fast (the prompt sound) and the out-of-step
        // motion rings on as the long tail, which is the treble's double decay (A4-F#5 early
        // 28 dB/s, late 3.7 against 33 and 3.2 recorded; at 1 cent the early slope was 16 and
        // the tail 6). No register is tuned closer or looser than another
        const double spread = unisonSpreadCents * (1.0 + 0.5 * keyNoise (k, 1));
        // Every key's pattern differs: which string is high or low, and by how much
        const double offsets3[3] = { -spread * (1.0 + 0.15 * keyNoise (k, 3)), 0.25 * spread * keyNoise (k, 2),
                                     spread * (1.0 + 0.15 * keyNoise (k, 4)) };
        const double flip = keyNoise (k, 5) < 0 ? -1.0 : 1.0;
        const double offsets2[2] = { -0.75 * flip * spread * (1.0 + 0.3 * keyNoise (k, 6)), 0.75 * flip * spread };
        // Centred, so the key's pitch stays where the tuner put it
        const double centre = key.numStrings == 3 ? (offsets3[0] + offsets3[1] + offsets3[2]) / 3.0
                            : key.numStrings == 2 ? (offsets2[0] + offsets2[1]) / 2.0 : 0.0;
        for (int s = 0; s < key.numStrings; ++s)
        {
            const double cents = (key.numStrings == 3 ? offsets3[s] : key.numStrings == 2 ? offsets2[s] : 0.0) - centre;
            const double f = key.targetFrequency * std::pow (2.0, cents / 1200.0);
            key.strings[(size_t) s] = makeString (designs[(size_t) k - 1], lengthOf (k), f);
        }

        const double n = k;
        // Hammers as measured on a Steinway D (Chabassier & Durufle, Inria RT-425, 2012): the
        // head's mass falls smoothly from 11 g at A0 to 5.8 g at C8, the felt gets harder
        // (stiffness K in F = K x^p) and its force law steeper (p from 2.3 to 3.5) up the
        // keyboard. No per-note voicing here: a technician's work is left to the voicing tables.
        key.hammerMass = -6.2348e-5 * n + 0.0112;
        key.feltExponent = 2.4295e-4 * n * n - 0.007703 * n + 2.337;
        // Chabassier's stiffness is used as measured. With a rigid hammer and no soft outer layer
        // (see ModalPiano::Settings) it gives the published hammer-string contact durations on a
        // D: about 2 ms at C4, 1.2 ms at C6 and 0.9 ms at C7 at mezzo-forte, 30% shorter at
        // fortissimo (Askenfelt & Jansson 1990). The old factor of 10 was compensating for the
        // outer layer and the hammer's internal springs, which lengthened every contact by 25-50%
        // and linearised the blow, so the tone brightened too little with velocity.
        key.feltStiffness = std::pow (10.0, 5.3097e-2 * n + 7.6425) * std::pow (1.0e-3, key.feltExponent);   // K per mm
        if (key.midiNote >= 84)
            key.feltStiffness *= kTopVoicing[std::min (24, key.midiNote - 84)];
        // Felt is not elastic: part of the push depends on how fast it is being squeezed, and
        // that part is lost as heat (Stulov 1995). Chabassier's D (RT-425, eq. 5, table A) gives
        // the ratio of the viscous to the elastic term as a time: 3.4 ms at C1, falling to
        // 0.3 ms at C3, 27 us at C5 and under 1 us at C8. So the bass felt loses much of each
        // blow (the contact lasts 3-4 ms there) and the treble felt is nearly elastic.
        key.feltDamping = std::pow (10.0, -4.366e-2 * n - 2.294);
        key.outerFeltDepth = lerp (1.0e-3, 0.1e-3, std::pow ((k - 1) / 87.0, 0.7));
        key.strikeRatio = kStrikeDistance[k - 1] / kLength[k - 1];   // about 0.12 of the speaking length on the D, top to bottom
        key.hammerWidth = lerp (12e-3, 5e-3, (k - 1) / 87.0);
        // No two hammers are alike: felt density, shape and wear differ a little from note
        // to note even on a well-voiced piano (which is what a voicer evens out by ear)
        key.hammerMass *= 1.0 + 0.03 * keyNoise (k, 11);
        key.feltStiffness *= 1.0 + 0.12 * keyNoise (k, 12);
        key.feltExponent += 0.08 * keyNoise (k, 13);
        key.outerFeltDepth *= 1.0 + 0.15 * keyNoise (k, 14);
        key.strikeRatio *= 1.0 + 0.025 * keyNoise (k, 15);
        key.hammerWidth *= 1.0 + 0.1 * keyNoise (k, 16);
        key.hasDamper = key.midiNote <= kLastDamperMidi;
        key.bassBridge = k <= kLastBassBridgeKey;
    }

    return keys;
}

} // namespace openpiano::piano
