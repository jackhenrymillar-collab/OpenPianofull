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
// Speaking lengths: Yamaha C5 (about 200 cm), the piano of the Salamander reference recordings.
//   No full C5 table is published, so these are the Steinway D lengths of Chabassier & Duruflé
//   (INRIA RT-0425, 2012) kept from C4 up and shortened below it: x0.80 at F2 rising to x1 at C4,
//   x0.66-0.72 on the bass bridge (A0 about 1.33 m), in line with measured 2 m grands.
//   The jump at key 21 is the break between the bass bridge and the treble bridge.
// String layout: monochords A0-E1, bichords F1-E2 (bass bridge), trichords F2-C8, with the
//   first trichords still copper-wound. Dampers stop at F6, as on concert grands.
// Hammers (key number n = 1..88), from J. O. Smith, "Physical Audio Signal Processing",
//   after Conklin, Hall and Stulov:
//   mass m = 11.074 - 0.074 n + 0.0001 n^2 grams
//   felt force F = Q0 * x^p, x = felt compression in mm, Q0 = 183 e^(0.045 n), p = 2.5 to 4 from bass to treble
// Wire: steel E = 2.0e11 Pa, density 7850 kg/m^3; copper winding 8960 kg/m^3.
//   Plain strings use standard piano-wire gauges in 0.025 mm steps, chosen for a tension of
//   690 N falling to 560 N at the top, fitted to the stiffness of a recorded grand
//   (about 1.0 mm at C4 and 0.8 mm at C8). Inharmonicity B = pi^3 E d^4 / (64 T L^2) from the
//   steel core, the standard result for stiff strings (Fletcher 1964).
// Tuning: A4 = 440 Hz. A3-A4 is the temperament octave; every other octave is tuned
//   upwards and downwards as a beatless octave on the actual (inharmonic) partials, as an
//   aural tuner would: 4:2 in the middle, 6:3 and 8:4 in the bass, 2:1 in the top octaves.
//   This produces the stretched "Railsback" curve of a real piano.

namespace
{
    constexpr double kSteinwayD[88] = {
        2.010, 2.009, 2.008, 2.007, 1.997, 1.981, 1.965, 1.938, 1.911, 1.879, 1.842, 1.805,
        1.762, 1.709, 1.655, 1.602, 1.548, 1.495, 1.442, 1.378, 1.837, 1.757,
        1.660, 1.591, 1.482, 1.403, 1.329, 1.259, 1.192, 1.129, 1.070, 1.013, 0.960, 0.909,
        0.861, 0.816, 0.773, 0.732, 0.694, 0.657, 0.622, 0.590, 0.559, 0.529,
        0.501, 0.475, 0.450, 0.426, 0.404, 0.383, 0.363, 0.344, 0.326, 0.308, 0.292, 0.277,
        0.262, 0.249, 0.236, 0.223, 0.211, 0.200, 0.190, 0.180, 0.171, 0.162,
        0.153, 0.145, 0.138, 0.130, 0.124, 0.117, 0.111, 0.105, 0.100, 0.095, 0.090, 0.085,
        0.081, 0.076, 0.072, 0.069, 0.065, 0.062, 0.058, 0.055, 0.052, 0.049
    };

    // Scale from Steinway D to Yamaha C5 lengths (see the source note above)
    double c5Scale (int key)
    {
        if (key >= 40) return 1.0;
        if (key > 20) return 0.80 + 0.20 * (key - 21) / 19.0;
        return 0.66 + 0.06 * (key - 1) / 19.0;
    }

    double lengthOf (int key) { return kSteinwayD[key - 1] * c5Scale (key); }

    constexpr double kSteelE = 2.0e11;
    constexpr double kSteelDensity = 7850.0;
    constexpr double kCopperDensity = 8960.0;
    constexpr double kWindingPacking = 0.85;

    constexpr int kLastMonochordKey = 8;   // E1
    constexpr int kLastBassBridgeKey = 20; // E2
    constexpr int kLastWoundKey = 26;      // A#2
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
        Design d {};
        if (key <= kLastWoundKey)
        {
            d.wound = true;
            // Core gauges, fitted so the inharmonicity follows a recorded grand (fitted B from
            // Salamander C5 single notes): stiffest at A0, falling to a minimum around F#2-A2
            // and rising again into the plain strings. Thick cores at the very bottom are what
            // give the lowest notes their slightly "out of tune with themselves" partials.
            constexpr double coreAt[] = { 2.03, 1.86, 1.72, 1.51, 1.45, 1.37, 1.22,  // keys 1, 4, ... 19
                                          1.27, 1.16, 1.12 };                      // keys 22, 25, 26 (treble bridge)
            if (key <= kLastBassBridgeKey)
            {
                const double t = (key - 1) / double (kLastBassBridgeKey - 1);
                const double pos = std::min ((key - 1) / 3.0, 5.999);
                const int i = (int) pos;
                d.coreDiameter = lerp (coreAt[i], coreAt[i + 1], pos - i) * 1e-3;
                d.woundTension = lerp (1800.0, 1350.0, t);
            }
            else
            {
                const double t = (key - kLastBassBridgeKey - 1) / double (kLastWoundKey - kLastBassBridgeKey - 1);
                const double pos = std::min ((key - 22) / 3.0, 1.0);
                d.coreDiameter = (key < 22 ? coreAt[7] : key >= 26 ? coreAt[9] : lerp (coreAt[7], coreAt[8], pos)) * 1e-3;
                d.woundTension = lerp (1150.0, 1000.0, t);
            }
            // The gauges were fitted on Steinway D lengths. B ~ d^4 / L^2, so keep the fitted
            // inharmonicity on the shorter C5 strings with a core thinner by sqrt(scale).
            d.coreDiameter *= std::sqrt (c5Scale (key));
            return d;
        }

        // Plain wire: pick the standard gauge that gives the target tension
        d.wound = false;
        const double t = (key - kLastWoundKey - 1) / double (88 - kLastWoundKey - 1);
        // Tension falls from 690 N to 560 N at the top: with the gauges this gives, the plain strings'
        // stiffness matches the recorded grand within about 10% (fitted with tools/calibrate/strings.py)
        d.plainTargetTension = lerp (690.0, 560.0, t);
        const double mu = d.plainTargetTension / std::pow (2.0 * length * frequency, 2);
        double diameter = std::sqrt (4.0 * mu / (M_PI * kSteelDensity));
        diameter = std::round (diameter / 0.025e-3) * 0.025e-3;     // standard gauges
        diameter = std::clamp (diameter, 0.775e-3, 1.25e-3);
        if (previousPlainDiameter > 0)
            diameter = std::min (diameter, previousPlainDiameter);  // gauges only get thinner going up
        d.coreDiameter = diameter;
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
        // Unisons drift apart after tuning: like the recorded Yamaha, a key's strings end up 1 to 3
        // cents apart (each outer string 0.5 to 1.5 cents from the centre; the two strings of a
        // bass bichord a little closer, about 0.75 to 2.3 cents, as wider ones pull its pitch off)
        // From A4 up the unisons sit closer (half as far apart up to A#5, three quarters to B6,
        // 0.3 in the top octave, and a quarter from A5 to C#6, where the recorded grand beats slowest;
        // 0.7 from E2 to G#4, which keeps their shimmer but beats no faster than the recording).
        // Close strings share the bridge motion strongly: the in-step vibration drains fast and the
        // out-of-step one is left ringing, which gives the recorded notes there their fast first
        // drop and long tail
        const double close = k >= 61 && k <= 66 ? 0.25 : k >= 28 && k < 50 ? 0.45 : k >= 20 && k < 49 ? 0.7 : k >= 49 && k < 63 ? 0.5 : k >= 63 && k < 67 ? 0.75 : k >= 67 && k < 75 ? 0.3 : k >= 75 ? 0.3 : 1.0;
        double spread = unisonSpreadCents * (1.0 + 0.5 * keyNoise (k, 1)) * close * (k < 43 ? 0.6 : 1.0);
        // (below D#5 the strings sit 40% closer than before: with the body coupled at its full,
        // physical strength, wider unisons there drain the long tail twice as fast as the recording)
        // Treble unisons are tuned by ear and some are left a little loose: on the recorded grand
        // about one key in four from D#5 up beats 3 to 5 times faster than its neighbours
        // The rest are a little looser too (Jack's ear: the loose ones sound better)
        if (k >= 43 && keyNoise (k, 7) > 0.5)
            spread *= 3.0 + 2.0 * (0.5 + 0.5 * keyNoise (k, 8));
        else if (k >= 43)
            spread *= 1.8;
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
        key.hammerMass = (11.074 - 0.074 * n + 0.0001 * n * n) * 1e-3;
        key.feltStiffness = 183.0 * std::exp (0.045 * n);
        // Hammer voicing for loud playing: a technician hardens or needles each felt so the
        // fortissimo brightness is even; fitted note by note to the recorded grand (every 3 keys from A0)
        {
            static constexpr double hardness[30] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5,
                                                     0.5, 2.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.67, 0.5,
                                                     0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
            const double pos = std::clamp ((k - 1) / 3.0, 0.0, 28.999);
            const int i = (int) pos;
            // (a lighter, harder head leaves the string sooner: brighter)
            key.hammerMass /= hardness[i] + (hardness[i + 1] - hardness[i]) * (pos - i);
        }
        key.feltExponent = 2.5 + 1.5 * (n - 1) / 87.0; // softer-curving felt: 2.5 in the bass to 4 in the treble
        key.outerFeltDepth = lerp (1.0e-3, 0.1e-3, std::pow ((k - 1) / 87.0, 0.7));
        key.strikeRatio = k <= 40 ? lerp (0.122, 0.115, (k - 1) / 39.0) : lerp (0.115, 0.08, (k - 40) / 48.0);
        key.hammerWidth = lerp (12e-3, 5e-3, (k - 1) / 87.0);
        // No two hammers are alike: felt density, shape and wear differ a little from note
        // to note even on a well-voiced piano (which is what a voicer evens out by ear)
        key.hammerMass *= 1.0 + 0.03 * keyNoise (k, 11);
        // C5 hammers are smaller than concert-grand ones in the bass, to suit the shorter
        // strings: about 15% lighter at A0, the same as the table from C4 up (estimate).
        key.hammerMass *= lerp (0.85, 1.0, std::clamp ((k - 1) / 39.0, 0.0, 1.0));
        // From the tenor up the hammer acts lighter on the string (its felt tip and head move
        // first, the rest follows through the shank), so the contact is short enough for loud
        // treble notes to come out as bright as recorded ones (half the mass from B3 up)
        key.hammerMass *= 1.0 - 0.5 * std::clamp ((k - 27) / 12.0, 0.0, 1.0);
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
