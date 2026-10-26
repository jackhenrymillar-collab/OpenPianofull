/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#include "ModalPiano.h"

#include <algorithm>
#include <cmath>
#include <complex>

namespace openpiano::piano
{

namespace
{
#include "BoardModes.inc"
    using cd = std::complex<double>;

    constexpr double kSteelE = 2.0e11;

    // Deterministic pseudo-random numbers in [0, 1)
    struct Rng
    {
        uint32_t s;
        explicit Rng (uint32_t seed) : s (seed * 747796405u + 2891336453u) { next(); }
        double next()
        {
            s ^= s << 13; s ^= s >> 17; s ^= s << 5;
            return (s & 0xffffff) / double (0x1000000);
        }
        double symmetric() { return 2.0 * next() - 1.0; } // [-1, 1)
    };

    // Internal (string material and air) losses, as a decay rate in 1/s for a partial
    // at frequency f, on a string whose fundamental is f1. Lower strings ring longer, and
    // high partials die faster. Bridge losses come on top, from the coupling itself.
    double internalDecay (double f, double f1, bool wound)
    {
        const double octavesAboveA0 = std::log2 (std::max (f1, 27.5) / 27.5);
        // The top two octaves also lose more in the string and its terminations: recorded
        // notes there fall 45 dB within one or two seconds
        const double b1 = 0.157 + 0.367 * octavesAboveA0 / 7.25 + 0.12 * std::pow (std::max (0.0, octavesAboveA0 - 5.0), 2.0);
        // Air and internal friction take the upper partials much faster (Bank 2000; Chaigne & Askenfelt 1994)
        const double b2 = wound ? 2.4e-7 : 1.2e-7;
        // Above about 2 kHz the losses climb faster still. Measured on recorded notes, overtones
        // at 3 to 6 kHz die at about 20 to 40 dB/s, so this rise is gentle
        // (the top fundamentals are spared: their losses are mostly at the bridge)
        const double overtone = std::clamp (f / f1 - 1.0, 0.0, 1.0);
        // From the middle up the strings themselves lose a quarter less, which gives the long
        // quiet tail of recorded notes there (about 3 dB/s from C4 to B5)
        const double b1s = b1 * (1.0 - 0.25 * std::clamp (octavesAboveA0 - 2.5, 0.0, 1.0));
        return b1s + b2 * f * f + 0.11 * overtone * std::pow (f / 2000.0, 3);
    }

    // Extra decay with the felt damper resting on the string
    double damperDecay (double f, double strength)
    {
        return strength * 10.0 * (1.0 + f / 2000.0);
    }

    // One sample of a run of string modes driven by the bridge acceleration a.
    // Returns their total force on the bridge.
    float stepModes (float* __restrict zr, float* __restrict zi, const float* __restrict pr,
                     const float* __restrict pi, const float* __restrict bgr, const float* __restrict bgi,
                     const float* __restrict fk, int n, float a)
    {
        float f = 0;
        for (int j = 0; j < n; ++j)
        {
            const float r = zr[j], im = zi[j];
            const float ni = pr[j] * im + pi[j] * r + bgi[j] * a;
            zr[j] = pr[j] * r - pi[j] * im + bgr[j] * a;
            zi[j] = ni;
            f += fk[j] * ni;
        }
        return f;
    }

    // Sum of u[k] * u[k + lag] over the valid range
    float lagProduct (const float* __restrict u, int n, int lag)
    {
        float s = 0;
        for (int k = 0; k + lag < n; ++k)
            s += u[k] * u[k + lag];
        return s;
    }

}

void ModalPiano::Modes::clear()
{
    for (auto* v : { &zr, &zi, &pr, &pi, &prFree, &piFree, &prDamp, &piDamp, &gr, &gi, &bgr, &bgi, &fk, &beta,
                     &phi, &pt, &wa, &prOs, &piOs, &gOsr, &gOsi, &sigmaInternal, &omega, &alt })
        v->clear();
}

ModalPiano::ModalPiano() = default;

void ModalPiano::prepare (double newSampleRate, const Settings& newSettings)
{
    sampleRate = newSampleRate;
    Ts = 1.0 / sampleRate;
    settings = newSettings;
    keys = buildConcertGrandScale (settings.unisonSpreadCents);

    // Age. An old piano drifts out of tune, key by key and string by string within each
    // unison; its strings lose energy faster (corrosion, fatigue, worn bearings), its
    // overtones scatter slightly off their ideal places, its hammers are packed hard and
    // grooved, its soundboard has lost much of its crown and some stiffness, and its
    // dampers and action are worn (those two are set by PianoInstrument).
    ageFactor = std::clamp (settings.age, 0.0, 1.0);
    if (ageFactor > 0)
    {
        Rng rng (4242);
        for (auto& key : keys)
        {
            // Untuned for years, the stretch partly relaxes (the treble sags and the bass
            // creeps back up, by up to a quarter), the whole piano settles a little flat, and
            // neighbouring sections of the scale drift together as the frame and board move
            const double position = (key.midiNote - kLowestKey) / double (kHighestKey - kLowestKey);
            const double stretchDrift = -0.25 * key.stretchCents - 5.0 * position
                                        + 4.0 * std::sin (2.0 * M_PI * position * 3.6 + 1.3);
            const double keyDrift = ageFactor * (stretchDrift + 7.0 * rng.symmetric());
            for (int i = 0; i < key.numStrings; ++i)
            {
                const double cents = keyDrift + 3.0 * ageFactor * rng.symmetric();
                key.strings[(size_t) i].fundamental *= std::pow (2.0, cents / 1200.0);
            }
        }
    }

    // Bridge zones: three along the short bass bridge, thirteen along the long treble bridge
    keyZone.assign (keys.size(), 0);
    {
        int firstTreble = 0;
        while (firstTreble < (int) keys.size() && keys[(size_t) firstTreble].bassBridge)
            ++firstTreble;
        for (int k = 0; k < (int) keys.size(); ++k)
            keyZone[(size_t) k] = k < firstTreble ? (kBassZones * k) / std::max (1, firstTreble)
                                                   : kBassZones + ((kZones - kBassZones) * (k - firstTreble)) / std::max (1, (int) keys.size() - firstTreble);
    }
    for (int z = 0; z < kZones; ++z)
    {
        const double p = std::clamp (layoutPosition (z), 0.0, 7.0);
        const int i = std::min (6, (int) p);
        zoneCoupling[(size_t) z] = zoneCouplingFit[(size_t) i] + (zoneCouplingFit[(size_t) i + 1] - zoneCouplingFit[(size_t) i]) * (p - i);
        zoneResistance[(size_t) z] = zoneResistanceFit[(size_t) i] + (zoneResistanceFit[(size_t) i + 1] - zoneResistanceFit[(size_t) i]) * (p - i);
    }

    buildModes();
    for (auto& st : strings)
        st.zone = keyZone[(size_t) st.key];
    buildSoundboard();
    buildFrame();
    buildAction();

    // Output filters: gentle top-end roll-off on the direct bridge path, and a DC blocker
    lpCoeff = std::exp (-2.0 * M_PI * 3500.0 / sampleRate);
    dcCoeff = std::exp (-2.0 * M_PI * 20.0 / sampleRate);

    fullList.reserve (strings.size());
    runs.reserve (2 * strings.size());
    segmentScratch.reserve (2 * strings.size());
    hammerList.reserve (keys.size());
    longList.reserve (strings.size());
    reset();
}

int ModalPiano::addString (const StringSpec& spec, int keyIndex, Kind kind, double strikeRatio, double hammerWidth,
                           int maxModes, double maxFrequency, uint32_t seed)
{
    const auto& key = keys[(size_t) keyIndex];
    const int os = std::max (1, settings.oversampling);
    Rng rng (seed);

    // No two strings are quite alike: small differences in where the hammer hits, in
    // stiffness and in losses keep the sound from being too regular.
    // Hammers are never set in a perfectly straight line, so each key's hammer also meets
    // its strings a little off the design point (shared by the unison)
    Rng keyRng ((uint32_t) keyIndex * 2654435761u + 77u);
    const double strike = strikeRatio * (1.0 + 0.025 * keyRng.symmetric()) * (1.0 + 0.012 * rng.symmetric());
    // Strings of one unison come off the same coil, but wire gauge varies by a fraction of a
    // percent, and stiffness goes with its fourth power: B differs by a few percent. That keeps
    // the upper partials of a unison from all beating in step, which would dull the tone at once.
    const double B = spec.inharmonicity * (1.0 + 0.04 * rng.symmetric());
    const double lossScale = (1.0 + 0.1 * rng.symmetric()) * (1.0 + 1.2 * ageFactor);
    const double fFlex = spec.fundamental / std::sqrt (1.0 + B);
    const double L = spec.length;

    String s;
    s.key = keyIndex;
    s.kind = kind;
    s.first = (int) m.size();
    s.strikeBase = (float) strike;
    s.halfWidthOverL = (float) (hammerWidth / (2.0 * L));
    s.mass = (float) (spec.linearDensity * L / 2.0);

    // Pan with the keyboard, as heard from the player's seat
    const double position = (key.midiNote - kLowestKey) / double (kHighestKey - kLowestKey);
    const double angle = M_PI / 4 + (position - 0.5) * 0.9 + 0.03 * rng.symmetric();
    // The direct (resistive) radiation makes up for the weaker board coupling at the top,
    // so the treble rings longer without getting quieter
    const double boost = 1.0 / zoneCoupling[(size_t) keyZone[(size_t) keyIndex]];
    s.panLeft = (float) (std::cos (angle) * boost);
    s.panRight = (float) (std::sin (angle) * boost);

    int sympathetic = 0;
    for (int n = 1; n <= maxModes; ++n)
    {
        const double f = n * fFlex * std::sqrt (1.0 + B * n * n)
                         * (n > 1 ? std::pow (2.0, 2.0 * ageFactor * rng.symmetric() / 1200.0) : 1.0);
        if (f > maxFrequency && n > 1)
            break;
        const double w = 2.0 * M_PI * f;
        const double sigma = internalDecay (f, spec.fundamental, spec.wound) * lossScale;
        const double sign = (n % 2 == 1) ? 1.0 : -1.0; // (-1)^(n+1)

        const cd pole = std::exp (cd (-sigma, w) * Ts);
        const cd poleOs = std::exp (cd (-sigma, w) * (Ts / os));
        const cd gain = pole * (Ts / w);          // impulse-invariant acceleration input
        const cd gainOs = poleOs * (Ts / os / w);

        // Bridge coupling. With the bridge end moving by u_b, mode n feels an
        // acceleration beta * u_b'' (beta = -2 (-1)^(n+1) / (n pi)), and the string pushes
        // on the bridge with T * dy/dx = T * (-1)^(n+1) * n pi / L per unit of mode amplitude.
        const double beta = -2.0 * sign / (n * M_PI);

        // Hammer: mode shape at the striking point, averaged over the felt's width
        const double x = n * M_PI * hammerWidth / (2.0 * L);
        const double width = std::abs (x) < 1e-9 ? 1.0 : std::sin (x) / x;

        m.zr.push_back (0); m.zi.push_back (0);
        m.pr.push_back ((float) pole.real()); m.pi.push_back ((float) pole.imag());
        m.prFree.push_back ((float) pole.real()); m.piFree.push_back ((float) pole.imag());
        m.prDamp.push_back (0); m.piDamp.push_back (0);
        m.gr.push_back ((float) gain.real()); m.gi.push_back ((float) gain.imag());
        m.bgr.push_back ((float) (gain.real() * beta)); m.bgi.push_back ((float) (gain.imag() * beta));
        m.fk.push_back ((float) (spec.tension * sign * n * M_PI / L));
        m.beta.push_back ((float) beta);
        m.alt.push_back ((float) sign);
        m.phi.push_back ((float) (std::sin (n * M_PI * strike) * width));
        m.pt.push_back ((float) std::sin (n * M_PI * strike));
        m.wa.push_back ((float) (-x * x / 6.0));
        m.prOs.push_back ((float) poleOs.real()); m.piOs.push_back ((float) poleOs.imag());
        m.gOsr.push_back ((float) gainOs.real()); m.gOsi.push_back ((float) gainOs.imag());
        m.sigmaInternal.push_back ((float) sigma);
        m.omega.push_back ((float) w);

        ++s.count;
        if (f < 4000.0)
            ++sympathetic;
    }

    // While only resonating, a string keeps its strongest low partials
    // (Undamped top strings have only a few modes and keep them all.)
    const bool segment = kind == Kind::Duplex || kind == Kind::Front;
    const int sympatheticCap = segment || ! key.hasDamper ? s.count : kind == Kind::Horizontal ? 6 : 16;
    s.sympatheticCount = ! key.hasDamper || segment ? s.count
                                                                 : std::clamp (sympathetic, 1, std::min (s.count, sympatheticCap));

    // Voicing: a technician evens out a piano's loudness note by note (hammer needling,
    // string levelling). This curve was fitted to the recorded grand's level, note by note
    // relative to the middle of the range. Pushing harder on the bridge while feeling its
    // motion less by the same factor keeps the decay and coupling unchanged.
    {
        static constexpr double voicingDb[30] = { -1.8, -0.6, -3.9, -9.8, -3.7, -11.4, -17.2, -17.6, -10.4, -9.1,
                                                  -1.1, -12.0, -4.9, 6.3, 7.7, 6.4, 1.3, 5.5, 7.7, 16.8,
                                                  -3.5, 5.0, 10.6, 9.7, 10.2, 6.9, 0.4, 3.0, 2.1, 7.7 }; // every 3 keys from A0
        const double pos = std::clamp ((key.midiNote - kLowestKey) / 3.0, 0.0, 28.999);
        const int i = (int) pos;
        const double db = voicingDb[i] + (voicingDb[i + 1] - voicingDb[i]) * (pos - i);
        const float g = (float) std::pow (10.0, db / 20.0);
        // In the top octave and a half the bridge and board give way less and less to the upper
        // partials (up to 15 dB less per octave of partial number from C7 up). Without this the
        // 2nd and 3rd partials there sit 10 to 20 dB above the recorded grand's: thin and metallic.
        const double topRoll = 15.0 * std::clamp ((keyIndex - 66) / 9.0, 0.0, 1.0);
        for (int j = 0; j < s.count; ++j)
        {
            const size_t q = (size_t) (s.first + j);
            const float gj = g * (float) std::pow (10.0, -topRoll * std::log2 (j + 1.0) / 20.0);
            m.fk[q] *= gj; m.bgr[q] /= gj; m.bgi[q] /= gj;
        }
    }

    strings.push_back (s);
    return (int) strings.size() - 1;
}

void ModalPiano::addLongitudinal (String& s, const StringSpec& spec)
{
    // Lengthwise waves travel at c = sqrt(EA / mu): the steel core carries the stretch and
    // the copper winding only adds mass. Mode m sits at m c / 2L, typically 12 to 20 times
    // the note's fundamental. The stretch is driven by the square of the string's slope,
    // so it is quiet in soft playing and grows quickly with force (Bank & Sujbert 2005).
    const double L = spec.length;
    const double EA = kSteelE * M_PI * spec.coreDiameter * spec.coreDiameter / 4.0;
    const double c2 = EA / spec.linearDensity;
    const double f1 = std::sqrt (c2) / (2.0 * L);
    const double fTop = std::min (0.4 * sampleRate, 7000.0);

    int partials = 0;
    for (int j = 0; j < s.count && j < settings.longitudinalPartials; ++j)
        if (m.omega[(size_t) modeIndex (s, j)] < 2.0 * M_PI * 5000.0)
            partials = j + 1;

    s.longFirst = (int) lm.zr.size();
    s.longCount = 0;
    s.longPartials = partials;
    s.tensionCoef = (float) (EA / 4.0 * std::pow (M_PI / L, 2));
    s.restTension = (float) spec.tension;

    for (int order = 1; order <= settings.longitudinalOrders && order * f1 < fTop && partials > 1; ++order)
    {
        const double w = 2.0 * M_PI * order * f1;
        // Free longitudinal ringing dies fast (end losses at the bridge pin and agraffe); what
        // remains audible is mostly the forced, phantom part that follows the transverse motion
        // Higher lengthwise modes lose energy faster still (end friction at the agraffe and
        // bridge pins, and in the winding), or they ring on as a thin, metallic whistle
        const double sigma = settings.longitudinalLoss + settings.longitudinalLossSlope * order * f1;
        const cd pole = std::exp (cd (-sigma, w) * Ts);
        const cd gain = pole * (Ts / w);
        lm.zr.push_back (0); lm.zi.push_back (0);
        lm.pr.push_back ((float) pole.real()); lm.pi.push_back ((float) pole.imag());
        lm.gr.push_back ((float) gain.real()); lm.gi.push_back ((float) gain.imag());
        lm.drive.push_back ((float) (-c2 * (order * M_PI / L) * std::pow (M_PI / L, 2) / 4.0));
        lm.out.push_back ((float) (EA * (order % 2 == 0 ? 1.0 : -1.0) * order * M_PI / L));
        lm.invW2.push_back ((float) (1.0 / (w * w)));
        lm.order.push_back (order);
        ++s.longCount;
    }
    scratchU.resize (std::max (scratchU.size(), (size_t) settings.longitudinalPartials));
}

void ModalPiano::buildModes()
{
    strings.clear();
    m.clear();
    lm = {};
    keyState.assign (keys.size(), {});
    alwaysOn.clear();

    const double fMax = std::min (0.45 * sampleRate, 16000.0);

    // Vertical motion: the hammer strikes it and it drives the soundboard
    for (size_t k = 0; k < keys.size(); ++k)
    {
        const auto& key = keys[k];
        auto& ks = keyState[k];
        ks.firstString = (int) strings.size();
        ks.numStrings = key.numStrings;
        for (int i = 0; i < key.numStrings; ++i)
        {
            const auto& spec = key.strings[(size_t) i];
            const int index = addString (spec, (int) k, Kind::Vertical, key.strikeRatio, key.hammerWidth,
                                         spec.wound ? 110 : 70, fMax, (uint32_t) (k * 7 + (size_t) i + 1));
            // The strings of a unison stretch almost identically, so one of them stands in
            // for all (its force is scaled up)
            (void) index;
        }
    }

    // Horizontal motion: tuned a hair differently (the bridge is stiffer sideways), barely
    // touched by the hammer, and fed mostly through the bridge
    if (settings.horizontal)
    {
        for (size_t k = 0; k < keys.size(); ++k)
        {
            const auto& key = keys[k];
            auto& ks = keyState[k];
            ks.firstHorizontal = (int) strings.size();
            for (int i = 0; i < key.numStrings; ++i)
            {
                Rng rng ((uint32_t) (k * 31 + (size_t) i * 5 + 17));
                auto spec = key.strings[(size_t) i];
                spec.fundamental *= std::pow (2.0, 0.25 * rng.symmetric() / 1200.0);
                const int index = addString (spec, (int) k, Kind::Horizontal, key.strikeRatio, key.hammerWidth, 24,
                                             std::min (fMax, 6000.0), (uint32_t) (k * 7 + (size_t) i + 1));
                // Hammer twist: a few hammers sit visibly off-square and push the strings
                // sideways much more; wear makes it worse
                Rng twistRng ((uint32_t) k * 7919u + 3u);
                const double twist = std::pow (twistRng.next(), 3.0) * (0.12 + 0.08 * ageFactor);
                strings[(size_t) index].hammerShare = (float) (0.03 + 0.04 * rng.next() + twist);
                strings[(size_t) (ks.firstString + i)].horizontal = index;
            }
        }
    }

    // Duplex scale: the short lengths of string between the bridge and the hitch pins
    // (and in front of the agraffes) in the treble are left free to ring. They are tuned
    // roughly to harmonics of the note, so they add a shimmer to the sustain.
    if (settings.duplex)
    {
        for (size_t k = 0; k < keys.size(); ++k)
        {
            const auto& key = keys[k];
            if (key.numStrings < 3 || key.bassBridge || key.midiNote < 54)
                continue;
            Rng rng ((uint32_t) (k * 13 + 101));
            const double t = (key.midiNote - 54) / double (kHighestKey - 54);
            const double ratio = (3.0 - 1.0 * t) * std::pow (2.0, 15.0 * rng.symmetric() / 1200.0);
            auto spec = key.strings[1];
            spec.length /= ratio;
            spec.fundamental *= ratio;
            spec.inharmonicity *= ratio * ratio;
            const int index = addString (spec, (int) k, Kind::Duplex, 0.0, 0.0, 6, std::min (fMax, 10000.0),
                                         (uint32_t) (k * 7 + 5000));
            keyState[k].duplexString = index;
            // Real duplex lengths are partly damped by the cloth and the hitch-pin end
            forModes (strings[(size_t) index], strings[(size_t) index].count, [&] (size_t j)
            {
                const float extra = (float) std::exp (-4.0 * Ts);
                m.pr[j] *= extra; m.pi[j] *= extra; m.prFree[j] *= extra; m.piFree[j] *= extra;
                m.bgr[j] *= 0.6f; m.bgi[j] *= 0.6f;
                m.sigmaInternal[j] += 4.0f;
            });
            alwaysOn.push_back (index);
        }
    }

    // More free string lengths. Behind the bridge, every treble string continues to its hitch
    // pin; above about C5 those rear lengths ring (the duplex above), below they are mostly
    // muted with felt but still sound faintly. In front, between the agraffes or the pressure
    // bar and the tuning pins, a short length rings too; it is fed by the far end of the string
    // through the frame, not by the bridge.
    for (size_t k = 0; k < keys.size(); ++k)
    {
        const auto& key = keys[k];
        if (key.bassBridge)
            continue;
        Rng rng ((uint32_t) (k * 31 + 707));
        // Muted rear lengths in the middle of the compass
        if (settings.duplex && key.midiNote >= 40 && key.midiNote < 54)
        {
            auto spec = key.strings[0];
            const double ratio = 3.3 * std::pow (2.0, 20.0 * rng.symmetric() / 1200.0);
            spec.length /= ratio; spec.fundamental *= ratio; spec.inharmonicity *= ratio * ratio;
            const int index = addString (spec, (int) k, Kind::Duplex, 0.0, 0.0, 4, std::min (fMax, 8000.0), (uint32_t) (k * 7 + 6000));
            // felt-muted: they lose energy quickly
            const float extra = (float) std::exp (-25.0 * Ts);
            forModes (strings[(size_t) index], strings[(size_t) index].count, [&] (size_t j)
            {
                m.pr[j] *= extra; m.pi[j] *= extra; m.prFree[j] *= extra; m.piFree[j] *= extra;
                m.sigmaInternal[j] += 25.0f;
            });
            alwaysOn.push_back (index);
        }
        // Front lengths, from about C4 up
        if (settings.farEnd > 0 && key.midiNote >= 60)
        {
            auto spec = key.strings[0];
            const double t = (key.midiNote - 60) / double (kHighestKey - 60);
            const double ratio = (6.0 - 2.5 * t) * std::pow (2.0, 25.0 * rng.symmetric() / 1200.0);
            spec.length /= ratio; spec.fundamental *= ratio; spec.inharmonicity *= ratio * ratio;
            const int index = addString (spec, (int) k, Kind::Front, 0.0, 0.0, 4, std::min (fMax, 12000.0), (uint32_t) (k * 7 + 7000));
            auto& st = strings[(size_t) index];
            // Not on the bridge at all
            forModes (st, st.count, [&] (size_t j) { m.bgr[j] = 0; m.bgi[j] = 0; m.fk[j] = 0;
                const float extra = (float) std::exp (-6.0 * Ts);
                m.pr[j] *= extra; m.pi[j] *= extra; m.prFree[j] *= extra; m.piFree[j] *= extra;
                m.sigmaInternal[j] += 6.0f; });
            alwaysOn.push_back (index);
        }
    }

    relayoutModes();

    // The strings of a unison stretch almost identically, so one of them stands in for all
    // (its force is scaled up)
    if (settings.longitudinal)
        for (size_t k = 0; k < keys.size(); ++k)
            if (keys[k].midiNote <= 84)
                addLongitudinal (strings[(size_t) keyState[k].firstString], keys[k].strings[0]);

    setDamperStrength (damperStrength);
}

void ModalPiano::relayoutModes()
{
    // addString put each string's modes in one block; split them into a low region and a
    // high region (see modeIndex)
    std::vector<size_t> order;
    order.reserve (m.size());
    std::vector<int> lowStart (strings.size()), highStart (strings.size());
    for (size_t i = 0; i < strings.size(); ++i)
    {
        lowStart[i] = (int) order.size();
        for (int k = 0; k < strings[i].sympatheticCount; ++k)
            order.push_back ((size_t) (strings[i].first + k));
    }
    for (size_t i = 0; i < strings.size(); ++i)
    {
        highStart[i] = (int) order.size();
        for (int k = strings[i].sympatheticCount; k < strings[i].count; ++k)
            order.push_back ((size_t) (strings[i].first + k));
    }
    for (auto* v : { &m.zr, &m.zi, &m.pr, &m.pi, &m.prFree, &m.piFree, &m.prDamp, &m.piDamp, &m.gr, &m.gi, &m.bgr,
                     &m.bgi, &m.fk, &m.beta, &m.phi, &m.pt, &m.wa, &m.prOs, &m.piOs, &m.gOsr, &m.gOsi, &m.sigmaInternal, &m.omega, &m.alt })
    {
        std::vector<float> moved (order.size());
        for (size_t j = 0; j < order.size(); ++j)
            moved[j] = (*v)[order[j]];
        *v = std::move (moved);
    }
    for (size_t i = 0; i < strings.size(); ++i)
    {
        strings[i].first = lowStart[i];
        strings[i].highFirst = highStart[i];
    }
}

void ModalPiano::setDamperStrength (double scale)
{
    if (! m.omega.empty() && std::abs (scale - damperStrength) < 1e-3 && m.prDamp[0] != 0.0f)
        return;
    damperStrength = scale;

    for (size_t j = 0; j < m.omega.size(); ++j)
    {
        const double w = m.omega[j];
        const double sigma = m.sigmaInternal[j] + damperDecay (w / (2.0 * M_PI), damperStrength);
        const cd pole = std::exp (cd (-sigma, w) * Ts);
        m.prDamp[j] = (float) pole.real();
        m.piDamp[j] = (float) pole.imag();
    }

    for (auto& s : strings)
        if (s.damped)
            forModes (s, s.count, [&] (size_t j)
            {
                m.pr[j] = m.prDamp[j];
                m.pi[j] = m.piDamp[j];
            });
}

double ModalPiano::layoutPosition (int zone)
{
    return zone < kBassZones ? (zone + 0.5) * 2.0 / kBassZones - 0.5
                             : 2.0 + (zone - kBassZones + 0.5) * 6.0 / (kZones - kBassZones) - 0.5;
}

void ModalPiano::buildSoundboard()
{
    // The soundboard: modes of a grand's spruce board, computed by tools/soundboard/build_board.py
    // (BoardModes.inc). The board has a grand's outline (straight spine, curved bentside, narrow
    // tail) and is clamped all round in the rim; it tapers from about 9 mm in the middle to 6 mm
    // at the edge; ribs every 13 cm stiffen it along the strings, so above about 1 kHz its
    // vibration is trapped between them; and the two bridges are stiff beams glued on top. So
    // the modes are lopsided and crowd into some parts of the board, and each part of the
    // bridge feels them differently. Their shapes at the 16 bridge zones come from that model.
    // The crown and the strings' downbearing stiffen the board further (below). Each mode's
    // mass is set so that the average admittance across the bridge matches measured bridge
    // admittances (Wogram 1980; Conklin 1996; Giordano 1998).
    static_assert (kBoardZones == kZones, "regenerate BoardModes.inc for this many bridge zones");
    Rng rng (1234);
    constexpr int Z = kZones;
    struct Mode { double f, loss, peak, radiation; int m, n; std::array<double, Z> phi; };
    std::vector<Mode> modes;

    const double lossFactor = 0.02 * (1.0 + 0.8 * ageFactor);
    const double meanBoardAdmittance = 0.9e-4;
    const double fc = 48.0 * settings.crown * (1.0 - 0.7 * ageFactor);
    // Above 2 kHz the modes overlap so much that the board acts as a plain resistance (the
    // direct bridge path below), so only the modes up to there are simulated one by one
    const double fTop = std::min (2000.0, 0.45 * sampleRate);
    for (int i = 0; i < kBoardModeCount; ++i)
    {
        const auto& src = kBoardModes[i];
        const double flat = src.f;
        const double fm = std::sqrt (flat * flat + fc * fc);
        if (fm > fTop)
            continue;
        // Mode spacing around this one (the ribs crowd the modes above about 1 kHz)
        int lo = std::max (0, i - 4), hi = std::min (kBoardModeCount - 1, i + 4);
        const double spacing = 2.0 * M_PI * (kBoardModes[hi].f - kBoardModes[lo].f) / std::max (1, hi - lo);
        // Losses rise with frequency (ribs, glue joints and radiation), so the upper modes
        // overlap and neighbouring notes come out evenly instead of hitting sharp peaks
        const double modeLoss = 2.0 * lossFactor * (1.0 + 3.0 * fm / 1000.0);
        Mode md { fm, modeLoss, 0, 1.0, 1, i, {} };
        const double sigmaFlat = M_PI * flat * modeLoss;
        md.peak = std::min (0.9e-3, meanBoardAdmittance * spacing / (M_PI * sigmaFlat)) * flat / fm;
        for (int z = 0; z < Z; ++z)
            md.phi[(size_t) z] = src.phi[z];
        modes.push_back (md);
    }
    const size_t plateModes = modes.size();

    // The rim: a few heavy, lightly damped low resonances that move the whole board edge.
    // The inner rim holds the board's edge, so they reach the bridge most near its two ends
    // (the bass bridge and the top treble run close to the rim).
    for (double f : { 68.0, 104.0, 151.0, 207.0, 286.0, 410.0 })
    {
        Mode md { f * (1.0 + 0.05 * rng.symmetric()), 0.008, 0.9e-4, 0.5, 0, 0, {} };
        for (int z = 0; z < Z; ++z)
        {
            const double p = layoutPosition (z);
            const double edge = std::max (std::exp (-p / 1.5), std::exp (-(7.0 - p) / 1.5));
            md.phi[(size_t) z] = 0.6 + 0.5 * edge;
        }
        modes.push_back (md);
    }
    // The back braces: heavy wooden beams under the board, joined to the rim. Each bends at
    // its own few hundred Hz and moves the part of the board (and bridge) above it.
    {
        const double braceF[] = { 118.0, 176.0, 243.0, 338.0, 471.0, 655.0 };
        const double braceAt[] = { 0.8, 2.6, 4.3, 5.8, 1.7, 6.8 }; // centre, in zones
        for (int k = 0; k < 6; ++k)
        {
            Mode md { braceF[k] * (1.0 + 0.04 * rng.symmetric()), 0.012, 0.35e-4, 0.25, -2, k, {} };
            for (int z = 0; z < Z; ++z)
                md.phi[(size_t) z] = std::exp (-std::pow ((layoutPosition (z) - braceAt[k]) / 1.3, 2.0));
            modes.push_back (md);
        }
    }
    // The cast-iron plate: stiff and barely damped. It carries the strings' far ends (the
    // agraffes and hitch pins), so it rings faintly under everything; its bending shapes
    // move the bass and treble ends in opposite directions.
    {
        const double plateF[] = { 122.0, 187.0, 324.0, 517.0, 768.0, 1130.0 };
        for (int k = 0; k < 6; ++k)
        {
            Mode md { plateF[k] * (1.0 + 0.03 * rng.symmetric()), 0.002, 2.5e-5, 0.3, -1, k, {} };
            for (int z = 0; z < Z; ++z)
                md.phi[(size_t) z] = std::cos (M_PI * (k + 1) * (layoutPosition (z) + 0.5) / 8.0);
            modes.push_back (md);
        }
    }

    const size_t n = modes.size();
    for (auto* v : { &board.wr, &board.wi, &board.qr, &board.qi, &board.hr, &board.hi, &board.omega, &board.sigma,
                     &board.radLeft, &board.radRight, &board.volume })
        v->assign (n, 0.0f);
    board.phi.assign (n * Z, 0.0f);
    board.scratchF.assign (n, 0.0f);
    board.scratchV.assign (n, 0.0f);

    // Sound radiated towards two listening points: the far-field (Rayleigh) integral of each
    // board mode's shape (computed with the shapes, in BoardModes.inc). The two points sit above
    // the keyboard, one towards the bass, one towards the treble, as from the player's seat.
    std::vector<double> rawL (n), rawR (n), rawV (n);
    for (size_t i = 0; i < plateModes; ++i)
    {
        const auto& src = kBoardModes[modes[i].n];
        rawL[i] = src.radLeft;
        rawR[i] = src.radRight;
        rawV[i] = src.volume;
    }

    // Keep the overall tonal balance the earlier, tuned board had: within each octave, scale
    // the radiation so its average strength matches; the differences between modes and
    // between the ears (the directivity) are what the shapes add.
    auto normaliseBands = [&] (std::vector<double>& w, double target)
    {
        for (double lo = 20.0; lo < 20000.0; lo *= 2.0)
        {
            double sum = 0; int count = 0;
            for (size_t i = 0; i < n; ++i)
                if (modes[i].m > 0 && modes[i].f >= lo && modes[i].f < 2.0 * lo) { sum += w[i] * w[i]; ++count; }
            if (count == 0 || sum <= 0) continue;
            const double scale = target / std::sqrt (sum / count);
            for (size_t i = 0; i < n; ++i)
                if (modes[i].m > 0 && modes[i].f >= lo && modes[i].f < 2.0 * lo) w[i] *= scale;
        }
    };
    normaliseBands (rawL, 0.85);
    normaliseBands (rawR, 0.85);
    // The air under the board is driven by the net volume a mode pushes: strong for the
    // lowest piston-like modes, small for the rest
    {
        double sum = 0; int count = 0;
        for (size_t i = 0; i < n; ++i)
            if (modes[i].m > 0 && modes[i].f < 400.0) { sum += rawV[i] * rawV[i]; ++count; }
        const double scale = count > 0 && sum > 0 ? 0.85 / std::sqrt (sum / count) : 0.0;
        for (auto& v : rawV) v *= scale;
    }

    for (size_t i = 0; i < n; ++i)
    {
        const double f = modes[i].f;
        const double w = 2.0 * M_PI * f;
        const double sigma = M_PI * f * modes[i].loss;
        // Mass from the admittance this mode should add, averaged over the bridge
        double meanPhi2 = 0;
        for (int z = 0; z < Z; ++z)
            meanPhi2 += modes[i].phi[(size_t) z] * modes[i].phi[(size_t) z];
        meanPhi2 = std::max (0.05, meanPhi2 / Z);
        const double modalMass = meanPhi2 / (2.0 * sigma * modes[i].peak);

        const cd pole = std::exp (cd (-sigma, w) * Ts);
        const cd gain = pole * (Ts / (w * modalMass));
        board.omega[i] = (float) w;
        board.sigma[i] = (float) sigma;
        board.qr[i] = (float) pole.real();
        board.qi[i] = (float) pole.imag();
        board.hr[i] = (float) gain.real();
        board.hi[i] = (float) gain.imag();
        // The bridge's own mass (see buildBody) also shields the strings from board modes
        // above its corner, most in the heavy bass and tenor
        for (int z = 0; z < Z; ++z)
        {
            const double fc = 2000.0 * std::pow (1.22, layoutPosition (z)), ratio = modes[i].f / fc;
            board.phi[(size_t) z * n + i] = (float) (zoneCoupling[(size_t) z] * modes[i].phi[(size_t) z] / std::sqrt (1.0 + ratio * ratio));
        }

        if (modes[i].m > 0)
        {
            board.radLeft[i] = (float) rawL[i];
            board.radRight[i] = (float) rawR[i];
            board.volume[i] = (float) rawV[i];
        }
        else
        {
            // Rim and plate: small radiators, much the same in both ears
            const double l = (0.6 + 0.4 * rng.next()) * (rng.next() < 0.5 ? -1 : 1);
            board.radLeft[i] = (float) (l * modes[i].radiation);
            board.radRight[i] = (float) (l * (0.8 + 0.2 * rng.next()) * modes[i].radiation);
        }
    }

    // Overall level: a mezzo-forte middle C peaks around -20 dBFS, fortissimo around -10
    outputGain = 17.0;
    buildBody();
}

void ModalPiano::buildFrame()
{
    // Cast-iron plate as seen from the agraffes and the pressure bar: much stiffer than the
    // soundboard (admittance around 1e-5 m/s/N at its resonances), and lightly damped
    Rng rng (99);
    const double fs[] = { 96.0, 171.0, 263.0, 402.0, 611.0, 893.0, 1370.0, 2110.0 };
    const size_t n = 8;
    for (auto* v : { &frame.wr, &frame.wi, &frame.qr, &frame.qi, &frame.hr, &frame.hi, &frame.omega, &frame.sigma, &frame.rad })
        v->assign (n, 0.0f);
    for (size_t i = 0; i < n; ++i)
    {
        const double f = fs[i] * (1.0 + 0.03 * rng.symmetric());
        const double w = 2.0 * M_PI * f, sigma = M_PI * f * 0.006;
        const double peak = 1.5e-5;
        const double mass = 1.0 / (2.0 * sigma * peak);
        const cd pole = std::exp (cd (-sigma, w) * Ts);
        const cd gain = pole * (Ts / (w * mass));
        frame.omega[i] = (float) w; frame.sigma[i] = (float) sigma;
        frame.qr[i] = (float) pole.real(); frame.qi[i] = (float) pole.imag();
        frame.hr[i] = (float) gain.real(); frame.hi[i] = (float) gain.imag();
        // The plate itself radiates only a little, but it is a large, ringing metal surface
        frame.rad[i] = (float) ((0.5 + 0.5 * rng.next()) * (rng.next() < 0.5 ? -1 : 1) * 0.6);
    }
    frame.accel = frame.lastVel = 0;
    frame.knock = 0;
    frame.knockDecay = std::exp (-1.0 / (0.0004 * sampleRate));
}

void ModalPiano::runFrame (double& left, double& right)
{
    // Far-end force of every sounding string: T dy/dx at x = 0 is T n pi / L per unit of mode
    // n, the bridge-end coupling without its alternating sign. The frame's motion moves that
    // end back, which the modes feel as beta (-1)^(n+1) times its acceleration.
    const float a = (float) (frame.accel * settings.farEnd);
    double F = 0;
    for (int index : frameList)
    {
        const auto& st = strings[(size_t) index];
        // The frame's modes all lie below a few kHz: the lowest partials do the coupling
        const int count = std::min (8, st.sympatheticCount);
        const size_t j0 = (size_t) st.first;
        float* __restrict zr = m.zr.data() + j0;
        float* __restrict zi = m.zi.data() + j0;
        const float* __restrict fk = m.fk.data() + j0;
        const float* __restrict alt = m.alt.data() + j0;
        const float* __restrict beta = m.beta.data() + j0;
        const float* __restrict gr = m.gr.data() + j0;
        const float* __restrict gi = m.gi.data() + j0;
        float f = 0;
        for (int j = 0; j < count; ++j)
        {
            f += fk[j] * alt[j] * zi[j];
            const float g = beta[j] * alt[j] * a;
            zr[j] += gr[j] * g;
            zi[j] += gi[j] * g;
        }
        F += f;
    }
    F *= settings.farEnd;
    // The hammer's blow also travels through the strings' far ends into the iron plate: a
    // short, sharp push that sets the frame ringing for a moment (the metallic shock of a
    // hard treble note)
    F += frame.knock;
    frame.knock *= frame.knockDecay;

    double v = 3e-6 * F, out = 0;
    for (size_t i = 0; i < frame.omega.size(); ++i)
    {
        const float x = frame.wr[i], y = frame.wi[i];
        const float nr = frame.qr[i] * x - frame.qi[i] * y + frame.hr[i] * (float) F;
        const float ni = frame.qr[i] * y + frame.qi[i] * x + frame.hi[i] * (float) F;
        frame.wr[i] = nr; frame.wi[i] = ni;
        const float vel = frame.omega[i] * nr - frame.sigma[i] * ni;
        v += vel;
        out += frame.rad[i] * vel;
    }
    frame.accel = (v - frame.lastVel) * sampleRate;
    frame.lastVel = v;
    left += out;
    right += 0.9 * out;
}

void ModalPiano::buildBody()
{
    // Air under a grand's board is boxed in by the rim and the floor and leaks out at the
    // gap underneath: a few broad resonances between about 90 and 350 Hz that fill out the bass.
    airModes.clear();
    // The first is the whole volume breathing through the gap (a Helmholtz-like mode); the
    // rest are standing waves across the case, about 2.4 m long and 1.45 m wide, up to 700 Hz.
    // The higher ones are weaker: the open bottom lets them leak away and the board drives
    // them less evenly.
    auto addAir = [&] (double f, double q, double g)
    {
        const double w = 2.0 * M_PI * f / sampleRate;
        const double alpha = std::sin (w) / (2.0 * q);
        const double a0 = 1.0 + alpha;
        airModes.push_back ({ (float) (alpha / a0), (float) (-2.0 * std::cos (w) / a0), (float) ((1.0 - alpha) / a0),
                              0, 0, 0, 0, (float) g });
    };
    addAir (92.0, 4.0, 0.55);
    {
        Rng rng (777);
        const double lx = 2.4, ly = 1.45, c = 343.0;
        for (int p = 0; p <= 6; ++p)
            for (int q = 0; q <= 4; ++q)
            {
                if (p == 0 && q == 0)
                    continue;
                const double f = 0.5 * c * std::sqrt (std::pow (p / lx, 2) + std::pow (q / ly, 2));
                if (f > 700.0)
                    continue;
                // Odd-odd shapes push and pull at the same time across the board and couple weakly
                const double shape = ((p % 2) == 1 && (q % 2) == 1) ? 0.4 : 1.0;
                const double g = 0.42 * shape / (1.0 + f / 160.0) * (0.8 + 0.4 * rng.next());
                addAir (f * (1.0 + 0.03 * rng.symmetric()), 5.0 + f / 60.0, g);
            }
    }
    airPressure = 0;

    // Lid: a reflection off its underside reaches the listener 1 to 2 ms after the direct
    // sound. The lid is about 1.5 m long, so it only reflects well above a few hundred Hz,
    // and its felt-free lacquered wood loses a little of the very top.
    const int maxDelay = (int) (0.004 * sampleRate) + 2;
    lidBufLeft.assign ((size_t) maxDelay, 0.0f);
    lidBufRight.assign ((size_t) maxDelay, 0.0f);
    lidPos = 0;
    lidDelayLeft = (int) (0.0021 * sampleRate);
    lidDelayRight = (int) (0.0013 * sampleRate);
    lidHpCoeff = std::exp (-2.0 * M_PI * 400.0 / sampleRate);
    lidLpCoeff = std::exp (-2.0 * M_PI * 9000.0 / sampleRate);
    lidHpLeft = lidHpRight = lidHpInLeft = lidHpInRight = lidLpLeft = lidLpRight = 0;

    // The bridge is a wooden beam with its own mass. Above a frequency set by that mass and
    // the board's resistance, the strings see the bridge's inertia rather than the board, so
    // they couple less. The bass bridge and the tenor are heavy; towards the top the bridge
    // is narrower and lighter, so the corner rises.
    for (int z = 0; z < kZones; ++z)
    {
        const double fc = std::min (2000.0 * std::pow (1.22, layoutPosition (z)), 0.4 * sampleRate);
        bridgeMassCoeff[(size_t) z] = std::exp (-2.0 * M_PI * fc / sampleRate);
        bridgeMassState[(size_t) z] = 0;
    }

    auto peaking = [this] (double f, double gainDb, double q)
    {
        Biquad b;
        const double A = std::pow (10.0, gainDb / 40.0), w = 2.0 * M_PI * f / sampleRate;
        const double alpha = std::sin (w) / (2.0 * q), a0 = 1.0 + alpha / A;
        b.b0 = (float) ((1.0 + alpha * A) / a0); b.b1 = (float) (-2.0 * std::cos (w) / a0);
        b.b2 = (float) ((1.0 - alpha * A) / a0); b.a1 = b.b1; b.a2 = (float) ((1.0 - alpha / A) / a0);
        return b;
    };
    // Above the simulated modes the board is not perfectly smooth: overlapping groups of
    // modes still leave broad peaks and dips of a few dB, different in each direction
    highBoardLeft.clear(); highBoardRight.clear();
    {
        Rng rng (4242);
        const double fs[] = { 5000.0, 6100.0, 7500.0, 9200.0, 11300.0 };
        for (double f : fs)
        {
            if (f > 0.42 * sampleRate)
                break;
            const double g = (rng.next() < 0.5 ? -1.0 : 1.0) * (2.0 + 2.0 * rng.next());
            highBoardLeft.push_back (peaking (f * (1.0 + 0.05 * rng.symmetric()), g, 2.5));
            highBoardRight.push_back (peaking (f * (1.0 + 0.05 * rng.symmetric()), -0.6 * g + 1.5 * rng.symmetric(), 2.5));
        }
    }
    // A medium room (about 9 x 7 x 4 m) with the piano a little off centre and the player at
    // the keyboard. Each ear gets the floor bounce, the near and far side walls, the ceiling
    // and the back wall, then a few second-order paths; later and more often reflected paths
    // are quieter and duller (wood floor, plaster, curtains).
    {
        reflections.clear();
        const double msL[] = { 3.1, 7.4, 11.8, 15.2, 19.7, 24.3, 29.8, 36.1, 43.5 };
        const double msR[] = { 3.4, 8.9, 10.6, 16.7, 21.4, 25.9, 31.6, 38.2, 45.0 };
        int maxDelay = 1;
        for (int i = 0; i < 9; ++i)
            for (int ear = 0; ear < 2; ++ear)
            {
                const double ms = ear == 0 ? msL[i] : msR[i];
                const double gain = 0.42 * std::exp (-ms / 22.0) * (i % 3 == 1 ? -1.0 : 1.0);
                const double cutoff = 7000.0 / (1.0 + 0.35 * i);
                Reflection r;
                r.delay = std::max (1, (int) (ms * 1e-3 * sampleRate));
                r.fromRight = (i + ear) % 2;           // reflections cross over between the ears
                r.toRight = ear == 1;
                r.gain = (float) gain;
                r.lpCoeff = (float) std::exp (-2.0 * M_PI * cutoff / sampleRate);
                r.state = 0;
                reflections.push_back (r);
                maxDelay = std::max (maxDelay, r.delay);
            }
        roomBufLeft.assign ((size_t) maxDelay + 2, 0.0f);
        roomBufRight.assign ((size_t) maxDelay + 2, 0.0f);
        roomPos = 0;
    }

    // The lid is a large spruce panel that the sound hits; it rings faintly at its own modes
    lidPanel.clear(); lidPanelGain.clear();
    {
        const double fs[] = { 68.0, 112.0, 171.0, 247.0, 352.0, 488.0 };
        for (double f : fs)
        {
            const double w = 2.0 * M_PI * f / sampleRate, alpha = std::sin (w) / (2.0 * 18.0), a0 = 1.0 + alpha;
            Biquad b;
            b.b0 = (float) (alpha / a0); b.b1 = 0; b.b2 = (float) (-alpha / a0);
            b.a1 = (float) (-2.0 * std::cos (w) / a0); b.a2 = (float) ((1.0 - alpha) / a0);
            lidPanel.push_back (b);
            lidPanelGain.push_back (0.05f);
        }
    }
    // Sound from the board spreads evenly at low frequencies but beams above a couple of kHz,
    // mostly out past the open lid (the right ear here)
    beamCoeff = std::exp (-2.0 * M_PI * 2000.0 / sampleRate);
    beamLpLeft = beamLpRight = 0;

    const double msL[3] = { 0.53, 1.37, 2.71 }, msR[3] = { 0.61, 1.63, 3.07 };
    for (int i = 0; i < 3; ++i)
    {
        diffLeft[(size_t) i].buf.assign ((size_t) std::max (1, (int) (msL[i] * 1e-3 * sampleRate)), 0.0f);
        diffRight[(size_t) i].buf.assign ((size_t) std::max (1, (int) (msR[i] * 1e-3 * sampleRate)), 0.0f);
        diffLeft[(size_t) i].pos = diffRight[(size_t) i].pos = 0;
    }
}

void ModalPiano::buildAction()
{
    // Wooden action parts ring briefly at a few hundred Hz to a few kHz (Askenfelt 1993,
    // "Observations on the transient components of the piano tone"). Three kinds of hit:
    //   0 key bottom: the key lands on the keybed felt, a dull low thump
    //   1 hammer knock: shank and hammer flex at the strike, a brighter click
    //   2 damper return: the felt lands back on the strings, soft and muffled
    // The last three are the keybed, the legs and the floor: the key's landing travels down
    // through the keyframe into the case and legs, a low thump you feel as much as hear
    const double freqs[11] = { 95, 160, 270, 430, 690, 1100, 1750, 2800, 31, 47, 66 };
    const double q[11] = { 6, 7, 8, 9, 10, 10, 11, 12, 3, 4, 5 };
    //   3 damper release: the felt meets a ringing string and briefly buzzes against it, a
    //     soft mid-range "chuff" whose strength follows how much the string was moving
    const double shapes[4][11] = {
        { 1.0, 0.9, 0.7, 0.45, 0.25, 0.12, 0.05, 0.02, 1.1, 0.9, 0.7 },
        { 0.15, 0.25, 0.4, 0.6, 0.8, 0.9, 0.7, 0.5, 0.1, 0.1, 0.1 },
        { 0.4, 0.5, 0.45, 0.3, 0.15, 0.06, 0.02, 0.01, 0.2, 0.2, 0.15 },
        { 0.1, 0.3, 0.55, 0.7, 0.6, 0.4, 0.2, 0.08, 0.0, 0.0, 0.0 },
    };
    for (int i = 0; i < ActionBank::N; ++i)
    {
        const double w = 2.0 * M_PI * freqs[i];
        const cd pole = std::exp (cd (-w / (2.0 * q[i]), w) * Ts);
        action.pr[(size_t) i] = (float) pole.real();
        action.pi[(size_t) i] = (float) pole.imag();
        for (int t = 0; t < ActionBank::Types; ++t)
            action.weight[(size_t) t][(size_t) i] = (float) shapes[t][i];
    }
    actionHits.clear();
    actionHits.reserve (512);
}

void ModalPiano::queueAction (int type, int delaySamples, double amount, int keyIndex)
{
    if (settings.actionNoise <= 0 || actionHits.size() >= actionHits.capacity())
        return;
    const auto& s = strings[(size_t) keyState[(size_t) keyIndex].firstString];
    actionHits.push_back ({ std::max (0, delaySamples), type, (float) (amount * settings.actionNoise), s.panLeft, s.panRight });
}

void ModalPiano::runAction (double& force, double& left, double& right)
{
    for (size_t i = 0; i < actionHits.size();)
    {
        auto& h = actionHits[i];
        if (h.delay-- > 0) { ++i; continue; }
        for (size_t j = 0; j < (size_t) ActionBank::N; ++j)
            action.zr[j] += h.amount * action.weight[(size_t) h.type][j];
        // The bank radiates from roughly where the keys were hit
        action.panLeft += 0.5f * (h.panLeft - action.panLeft);
        action.panRight += 0.5f * (h.panRight - action.panRight);
        h = actionHits.back();
        actionHits.pop_back();
    }
    float f = 0;
    for (size_t j = 0; j < (size_t) ActionBank::N; ++j)
    {
        const float r = action.zr[j], im = action.zi[j];
        action.zr[j] = action.pr[j] * r - action.pi[j] * im;
        action.zi[j] = action.pr[j] * im + action.pi[j] * r;
        f += action.zi[j];
    }
    force += f;
    left += f * action.panLeft;
    right += f * action.panRight;
}

void ModalPiano::reset()
{
    for (auto& st : strings)
    {
        const float ratio = 1.0f / st.outScale;
        forModes (st, st.count, [&] (size_t q) { m.fk[q] *= ratio; m.bgr[q] /= ratio; m.bgi[q] /= ratio; });
        st.outScale = 1;
    }

    action.zr.fill (0); action.zi.fill (0);
    actionHits.clear();
    std::fill (m.zr.begin(), m.zr.end(), 0.0f);
    std::fill (m.zi.begin(), m.zi.end(), 0.0f);
    std::fill (lm.zr.begin(), lm.zr.end(), 0.0f);
    std::fill (lm.zi.begin(), lm.zi.end(), 0.0f);
    std::fill (board.wr.begin(), board.wr.end(), 0.0f);
    std::fill (board.wi.begin(), board.wi.end(), 0.0f);
    lastBridgeVelocity = bridgeAcceleration = 0;
    zoneForce.fill (0); zoneAcceleration.fill (0); lastZoneVelocity.fill (0); zoneBoardVelocity.fill (0);
    lastHorizontalVelocity = horizontalAcceleration = 0;
    resistiveLpLeft = resistiveLpRight = 0;
    buildBody();
    buildFrame();
    dcLeft = dcRight = dcInLeft = dcInRight = 0;
    sustain = false;
    sostenutoDown = false;

    for (auto& k : keyState)
    {
        k.down = false;
        k.sostenuto = false;
        k.hammer = {};
    }
    for (auto& s : strings)
    {
        s.activity = Activity::Off;
        s.damped = true;
    }
    // Strings without dampers (the top of the compass) listen to the bridge all the time,
    // whether or not their keys have been played: the sympathetic halo of the treble
    for (size_t k = 0; k < keys.size(); ++k)
        if (! keys[k].hasDamper)
            updateDamping ((int) k);
    // Duplex segments have no dampers and always listen
    for (int i : alwaysOn)
    {
        auto& s = strings[(size_t) i];
        s.damped = false;
        s.activity = Activity::Sympathetic;
        forModes (s, s.count, [&] (size_t j)
        {
            m.pr[j] = m.prFree[j];
            m.pi[j] = m.piFree[j];
        });
    }
    for (size_t k = 0; k < keys.size(); ++k)
        updateDamping ((int) k);
    checkActivity();
}

void ModalPiano::setActivity (String& s, Activity a)
{
    if (s.activity == a)
        return;

    // Modes that stop being computed are cleared, so they restart from rest
    const int keep = a == Activity::Full ? s.count : a == Activity::Sympathetic ? s.sympatheticCount : 0;
    int k = 0;
    forModes (s, s.count, [&] (size_t j)
    {
        if (k++ >= keep)
            m.zr[j] = m.zi[j] = 0.0f;
    });
    if (a != Activity::Full)
        for (int j = s.longFirst; j < s.longFirst + s.longCount; ++j)
            lm.zr[(size_t) j] = lm.zi[(size_t) j] = 0.0f;
    s.activity = a;
}

void ModalPiano::updateDamping (int keyIndex)
{
    const auto& key = keys[(size_t) keyIndex];
    auto& ks = keyState[(size_t) keyIndex];
    const bool damped = key.hasDamper && ! ks.down && ks.damperPending <= 0 && ! sustain && ! ks.sostenuto;

    forEachString (ks, [&] (String& s)
    {
        if (s.damped != damped || s.activity == Activity::Off)
        {
            s.damped = damped;
            const auto& prSrc = damped ? m.prDamp : m.prFree;
            const auto& piSrc = damped ? m.piDamp : m.piFree;
            forModes (s, s.count, [&] (size_t j)
            {
                m.pr[j] = prSrc[j];
                m.pi[j] = piSrc[j];
            });
        }
        // Undamped strings listen to the bridge
        if (! damped && s.activity == Activity::Off)
            setActivity (s, Activity::Sympathetic);
    });
}

void ModalPiano::setKeyDown (int midiNote, bool down)
{
    const int k = midiNote - kLowestKey;
    if (k < 0 || k >= (int) keys.size())
        return;
    auto& ks = keyState[(size_t) k];
    if (! down && ks.down)
    {
        // The key takes a moment to rise far enough for the damper to reach the strings.
        // Re-press it before then (a fast repeat or a half-stroke) and the damper never lands.
        ks.releasedAt = sampleClock;
        ks.damperPending = std::max (1, (int) (0.02 * sampleRate));
        ks.down = false;
        return;
    }
    ks.down = down;
    ks.damperPending = 0;
    updateDamping (k);
    checkActivity();
}

void ModalPiano::setSustain (bool down)
{
    if (sustain == down)
        return;
    sustain = down;
    for (size_t k = 0; k < keys.size(); ++k)
        updateDamping ((int) k);
    checkActivity();
}

void ModalPiano::setSostenuto (bool down)
{
    if (sostenutoDown == down)
        return;
    sostenutoDown = down;
    for (size_t k = 0; k < keys.size(); ++k)
    {
        auto& ks = keyState[k];
        // The sostenuto rod catches every damper that is lifted when the pedal goes down
        ks.sostenuto = down && keys[k].hasDamper && (ks.down || ks.damperPending > 0 || sustain);
        updateDamping ((int) k);
    }
    checkActivity();
}

void ModalPiano::strike (int midiNote, double hammerVelocity)
{
    const int k = midiNote - kLowestKey;
    if (k < 0 || k >= (int) keys.size())
        return;

    auto& ks = keyState[(size_t) k];
    // Repetition. A grand's double escapement lets a key be struck again before it has risen
    // all the way, but from a partly lifted key the hammer has less travel to gather speed.
    // Fully lifted after about 90 ms.
    if (! ks.down && ks.releasedAt >= 0)
    {
        const double lift = std::clamp ((sampleClock - ks.releasedAt) / (0.09 * sampleRate), 0.0, 1.0);
        hammerVelocity *= 0.55 + 0.45 * std::sqrt (lift);
    }
    ks.releasedAt = -1;
    ks.damperPending = 0;
    ks.down = true;
    updateDamping (k);
    forEachString (ks, [&] (String& s) { setActivity (s, Activity::Full); });

    // The hammer is never perfectly level, so it meets the strings of a unison a few
    // hundredths of a millimetre apart. This is what excites the slow "aftersound".
    Rng rng ((uint32_t) midiNote * 977u + 13u);
    double nearest = 1e9;
    for (int i = 0; i < ks.numStrings; ++i)
    {
        const auto& s = strings[(size_t) (ks.firstString + i)];
        ks.hammer.gap[(size_t) i] = settings.hammerUnevenness * rng.next();
        // The felt is not uniform across the hammer's width either: each string of the
        // unison is pushed a little harder or softer. This uneven start is what leaves the
        // strings swinging partly against each other, the slow, quiet "aftersound" tail.
        ks.hammer.feltScale[(size_t) i] = 1.0 + feltScatter * (2.0 * Rng ((uint32_t) (k * 3 + i) * 7919u + 11u).next() - 1.0);
        // Una corda: the whole action slides sideways. On a trichord the hammer misses the
        // first string, which is left to ring in sympathy; on a bichord it only grazes one
        // string. Every hammer also meets the strings with softer felt between its grooves.
        if (softPedal > 0)
        {
            const double miss = ks.numStrings == 3 ? 1.0 : ks.numStrings == 2 ? 0.6 : 0.0;
            if (i == 0 && miss > 0)
            {
                ks.hammer.feltScale[(size_t) i] *= 1.0 - miss * softPedal;
                ks.hammer.gap[(size_t) i] += miss * softPedal * 0.3e-3;
            }
        }
        double u = 0;
        forModes (s, s.count, [&] (size_t j) { u += m.phi[j] * m.zi[j]; });
        nearest = std::min (nearest, u + ks.hammer.gap[(size_t) i]);
    }

    // Action noises scale with how hard the key is played. The key lands on the keybed a
    // few milliseconds after the hammer reaches the strings; the knock is louder in the
    // treble, where hammers are smaller and the strings shorter.
    {
        const double v = std::max (0.05, hammerVelocity);
        const double treble = (keys[(size_t) k].midiNote - kLowestKey) / double (kHighestKey - kLowestKey);
        queueAction (0, (int) (0.003 * sampleRate), 0.4 * v * (1.0 - 0.5 * treble), k);
        // In the top octaves the knock is a large part of what you hear: a woody thud that
        // tells a piano from a harp or a glockenspiel
        const double top = std::pow (std::max (0.0, (treble - 0.6) / 0.4), 2.0);
        queueAction (1, 0, 0.09 * std::pow (v, 1.5) * (0.4 + 0.6 * treble + 1.2 * top), k);
        if (settings.actionNoise > 0)
            frame.knock += settings.actionNoise * settings.farEnd * frameKnockGain * std::pow (v, 1.5) * (0.2 + 0.8 * treble * treble);
    }

    // Soft notes are made a little louder without changing their tone: the key's strings push
    // harder on the bridge and feel its motion less by the same factor, so their decay and
    // coupling stay as they were. Fortissimo (5 m/s) is left unchanged.
    {
        // Recorded pianos have more dynamic range in the treble than in the bass, so the lift
        // fades out towards the top of the keyboard
        const double top = k / double (kHighestKey - kLowestKey);
        const float lift = (float) std::pow (5.0 / std::clamp (hammerVelocity, 0.3, 5.0), settings.softLift * (1.0 - top));
        forEachString (ks, [&] (String& s)
        {
            const float ratio = lift / s.outScale;
            forModes (s, s.count, [&] (size_t j) { m.fk[j] *= ratio; m.bgr[j] /= ratio; m.bgi[j] /= ratio; });
            s.outScale = lift;
        });
    }

    // Strike point moves a little with how hard the key is played: on a hard blow the shank
    // bends and the head rolls, so the felt lands slightly further from the agraffe
    {
        const double shift = 0.012 * (std::clamp (hammerVelocity, 0.3, 5.0) - 2.0) / 3.0;
        forEachString (ks, [&] (String& s)
        {
            if (s.kind != Kind::Vertical)
                return;
            const double strikePos = s.strikeBase * (1.0 + shift);
            forModes (s, s.count, [&] (size_t j)
            {
                const double n = std::round (2.0 / (M_PI * std::max (1e-6f, std::abs (m.beta[j]))));
                const double x = n * M_PI * s.halfWidthOverL;
                const double width = std::abs (x) < 1e-9 ? 1.0 : std::sin (x) / x;
                m.phi[j] = (float) (std::sin (n * M_PI * strikePos) * width);
                m.pt[j] = (float) std::sin (n * M_PI * strikePos);
            });
        });
    }

    // Let-off: the jack slips out from under the hammer a couple of millimetres before the
    // string, and the hammer flies the rest of the way free, slowed a little by gravity. Soft
    // notes therefore arrive a few milliseconds later than hard ones.
    const double letOff = 2.0e-3 * (1.0 + 0.15 * rng.symmetric()) * (1.0 + 0.5 * ageFactor);
    hammerVelocity = std::sqrt (std::max (0.0025, hammerVelocity * hammerVelocity - 2.0 * 9.81 * letOff));

    // Above about MIDI velocity 100 the felt is squeezed into its densest state and its
    // force law steepens further, so the hardest blows get disproportionately brighter
    ks.hammer.exponentBoost = 0.5 * std::clamp ((hammerVelocity - 3.3) / 1.7, 0.0, 1.0);
    ks.hammer.softFelt = softPedal;

    ks.hammer.active = true;
    ks.hammer.position = nearest - letOff; // just let off, flying free
    ks.hammer.velocity = std::max (0.05, hammerVelocity);
    ks.hammer.headPos = ks.hammer.shankPos = nearest - letOff;
    ks.hammer.headVel = ks.hammer.shankVel = ks.hammer.velocity;
    ks.hammer.samplesSinceStrike = 0;
    ks.hammer.samplesApart = 0;
    for (int i = 0; i < ks.numStrings; ++i)
    {
        const auto& s = strings[(size_t) (ks.firstString + i)];
        double u = 0;
        forModes (s, s.count, [&] (size_t j) { u += m.phi[j] * m.zi[j]; });
        ks.hammer.lastCompression[(size_t) i] = nearest - letOff - u - ks.hammer.gap[(size_t) i];
        ks.hammer.forceSum[(size_t) i] = 0;
        ks.hammer.feltHistory[(size_t) i] = 0;
    }
    checkActivity();
}

double ModalPiano::partialFrequency (int midiNote, int stringIndex, int k) const
{
    const int key = midiNote - kLowestKey;
    if (key < 0 || key >= (int) keys.size())
        return 0;
    const auto& ks = keyState[(size_t) key];
    if (stringIndex >= ks.numStrings)
        return 0;
    const auto& s = strings[(size_t) (ks.firstString + stringIndex)];
    if (k < 1 || k > s.count)
        return 0;
    return m.omega[(size_t) modeIndex (s, k - 1)] / (2.0 * M_PI);
}

void ModalPiano::checkActivity()
{
    // Retire strings that have gone quiet, and rebuild the work lists
    constexpr float silent = 1e-18f; // (1 nm)^2 displacement
    fullList.clear();
    hammerList.clear();
    longList.clear();
    frameList.clear();
    activeStrings = 0;

    auto visit = [&] (int index, bool hammerActive)
    {
        auto& s = strings[(size_t) index];
        if (s.activity == Activity::Off)
            return;

        if (! hammerActive)
        {
            // low: the modes kept while resonating; high: the rest (or all of them, if only
            // resonating, so a quiet resonating string can still be retired)
            const bool full = s.activity == Activity::Full;
            float low = 0, high = 0;
            int k = 0;
            forModes (s, full ? s.count : s.sympatheticCount, [&] (size_t j)
            {
                const float e = m.zr[j] * m.zr[j] + m.zi[j] * m.zi[j];
                if (full && k++ < s.sympatheticCount)
                    low = std::max (low, e);
                else
                    high = std::max (high, e);
            });

            if (s.damped && std::max (low, high) < silent)
                setActivity (s, Activity::Off);
            else if (! s.damped && s.activity == Activity::Full && high < silent && low < silent * 1e4f)
                setActivity (s, Activity::Sympathetic);
        }

        if (s.activity == Activity::Full && s.kind == Kind::Vertical)
            ++activeStrings;
        if (s.activity == Activity::Full && s.longCount > 0)
            longList.push_back (index);
        if (s.activity == Activity::Full && s.kind == Kind::Vertical)
            frameList.push_back (index);
        // Vertical strings under a hammer are stepped by the hammer code instead
        if (s.activity != Activity::Off && ! (hammerActive && s.kind == Kind::Vertical))
            fullList.push_back (index);
    };

    for (size_t k = 0; k < keys.size(); ++k)
    {
        const auto& ks = keyState[k];
        if (ks.hammer.active)
            hammerList.push_back ((int) k);
        for (int i = 0; i < ks.numStrings; ++i)
        {
            visit (ks.firstString + i, ks.hammer.active);
            if (ks.firstHorizontal >= 0)
                visit (ks.firstHorizontal + i, ks.hammer.active);
        }
    }
    for (int i : alwaysOn)
    {
        fullList.push_back (i);
        if (strings[(size_t) i].kind == Kind::Front)
            frameList.push_back (i);
    }

    // Merge the computed mode ranges into runs where they sit back to back, at most 12
    // strings (about 4 keys) each, so their pan barely differs
    auto& segments = segmentScratch;
    segments.clear();
    for (int index : fullList)
    {
        const auto& s = strings[(size_t) index];
        const bool side = s.kind == Kind::Horizontal;
        segments.push_back ({ s.first, s.sympatheticCount, s.panLeft, s.panRight, side, s.zone });
        if (s.activity == Activity::Full && s.count > s.sympatheticCount)
            segments.push_back ({ s.highFirst, s.count - s.sympatheticCount, s.panLeft, s.panRight, side, s.zone });
    }
    std::sort (segments.begin(), segments.end(), [] (const Segment& a, const Segment& b) { return a.first < b.first; });

    runs.clear();
    int runSegments = 0;
    for (const auto& g : segments)
    {
        if (! runs.empty() && runs.back().first + runs.back().count == g.first && runs.back().side == g.side
            && runs.back().zone == g.zone && runSegments < 12)
        {
            auto& r = runs.back();
            r.count += g.count;
            r.panLeft += (g.panLeft - r.panLeft) / float (runSegments + 1);
            r.panRight += (g.panRight - r.panRight) / float (runSegments + 1);
            ++runSegments;
        }
        else
        {
            runs.push_back ({ g.first, g.count, g.panLeft, g.panRight, g.side, g.zone });
            runSegments = 1;
        }
    }
}

void ModalPiano::runHammerKey (int keyIndex, double ub, double& forceSum, double& forceLeft, double& forceRight)
{
    auto& ks = keyState[(size_t) keyIndex];
    auto& h = ks.hammer;
    const auto& key = keys[(size_t) keyIndex];
    const int os = std::max (1, settings.oversampling);
    const double dt = Ts / os;
    // Felt grooves: years of playing pack the felt hard where the strings land, most in the
    // middle of the keyboard where most music sits. The grooves are harder and flatter than
    // the crown, so a worn hammer sounds brighter and more cutting there.
    const double groove = ageFactor * std::exp (-std::pow ((key.midiNote - 66) / 18.0, 2.0))
                          * (0.7 + 0.6 * Rng ((uint32_t) keyIndex * 131u + 9u).next());
    // The una corda moves the strings onto the unpacked felt beside the grooves: softer and
    // with a gentler force law, so the tone is darker as well as quieter
    const double q0 = key.feltStiffness * hardness * (1.0 + 1.5 * ageFactor) * (1.0 + 0.9 * groove * (1.0 - h.softFelt))
                      * (1.0 - 0.35 * h.softFelt);
    const double p = key.feltExponent + h.exponentBoost - 0.3 * h.softFelt;
    // Felt has memory (Stulov 1995): squeezed quickly it is stiff; squeezed slowly it has
    // time to relax and gives way. So a gentle touch makes a long, soft contact that only
    // excites the low partials, and a hard blow a short, stiff one that brings out the top.
    const double memory = settings.feltMemory;
    const double coreDepth = key.outerFeltDepth * 1000.0 * settings.outerFeltScale * (1.0 - 0.6 * ageFactor) * (1.0 - 0.5 * groove); // [mm]
    // Treble hammers are smaller and their outer felt firmer, or soft notes there would vanish
    const double outer = settings.outerFelt * std::pow (10.0, std::max (0.0, (keyIndex - 39) / 48.0) * trebleFirmness);
    const double relax = 1.0 - std::exp (-dt / std::max (1e-6, settings.feltRelaxTime));
    // Squeeze at which the contact patch is hammerWidth long (a firm mezzo-forte blow)
    const double contactRef = 1.55e-3 * std::exp (-0.0185 * keyIndex);

    bool touching = false;
    for (int step = 0; step < os; ++step)
    {
        double hammerForce = 0;
        for (int i = 0; i < ks.numStrings; ++i)
        {
            auto& s = strings[(size_t) (ks.firstString + i)];

            double u = 0;
            forModes (s, s.count, [&] (size_t j) { u += m.phi[j] * m.zi[j]; });

            const double compression = h.position - u - h.gap[(size_t) i];
            h.lastCompression[(size_t) i] = compression;
            // Two layers: a soft outer felt, then the hard, dense core once the outer is squeezed through
            const double xmm = compression * 1000.0;
            const double elastic = compression > 0
                ? outer * std::pow (xmm, p) + std::pow (std::max (0.0, xmm - coreDepth), p)
                : 0.0;
            auto& history = h.feltHistory[(size_t) i];
            history += (elastic - history) * relax;
            const double f = compression > 0 ? std::max (0.0, q0 * h.feltScale[(size_t) i] * (elastic - memory * history)) : 0.0;

            // The felt crown is round: squeezed by d, it flattens into a patch about
            // 2 sqrt(2 R d) long. A light touch meets the string almost at a point; a hard
            // blow spreads over a wider patch, which averages out the shortest wavelengths.
            // The patch is hammerWidth at the reference squeeze; the mode shape is
            // sin(n pi a) sinc(x), with sinc(x) = 1 - x^2/6 + x^4/120 (x stays below about 2).
            if (settings.contactSpread && compression > 0)
            {
                const float q = (float) std::clamp (compression / contactRef, 0.09, 2.56); // (width / reference)^2
                forModes (s, s.count, [&] (size_t j)
                {
                    const float a = m.wa[j] * q;
                    m.phi[j] = m.pt[j] * (1.0f + a + 0.3f * a * a);
                });
            }
            hammerForce += f;
            h.forceSum[(size_t) i] += f;
            touching = touching || f > 0;

            const float hammerAcc = (float) (f / s.mass);
            const float bridgeAcc = (float) ub / s.outScale;
            forModes (s, s.count, [&] (size_t j)
            {
                const float acc = m.beta[j] * bridgeAcc + m.phi[j] * hammerAcc;
                const float zr = m.zr[j], zi = m.zi[j];
                m.zr[j] = m.prOs[j] * zr - m.piOs[j] * zi + m.gOsr[j] * acc;
                m.zi[j] = m.prOs[j] * zi + m.piOs[j] * zr + m.gOsi[j] * acc;
            });
        }

        // Hammer flight (symplectic Euler). The hammer is not rigid: its felt tip sits on a
        // wooden head that rocks on the shank (about 1.5 kHz), and the shank bends (about
        // 400 Hz). During the first instants of contact only the tip and head push on the
        // string; the rest of the mass arrives through the flexing shank a little later,
        // which rounds off the force pulse (Askenfelt & Jansson 1991; Giordano & Millis 2004).
        if (settings.shankFlex > 0)
        {
            const double mt = 0.3 * key.hammerMass, mh = 0.45 * key.hammerMass, ms = 0.25 * key.hammerMass;
            const double fr = 1500.0 / std::sqrt (settings.shankFlex), fb = 400.0 / std::sqrt (settings.shankFlex);
            const double mr = mt * mh / (mt + mh), mb = mh * ms / (mh + ms);
            const double kr = std::pow (2.0 * M_PI * fr, 2) * mr, kb = std::pow (2.0 * M_PI * fb, 2) * mb;
            const double cr = 2.0 * 0.4 * std::sqrt (kr * mr), cb = 2.0 * 0.1 * std::sqrt (kb * mb);
            const double fRock = kr * (h.headPos - h.position) + cr * (h.headVel - h.velocity);
            const double fBend = kb * (h.shankPos - h.headPos) + cb * (h.shankVel - h.headVel);
            h.velocity += (fRock - hammerForce) / mt * dt;
            h.headVel += (fBend - fRock) / mh * dt;
            h.shankVel += (-fBend) / ms * dt;
            h.position += h.velocity * dt;
            h.headPos += h.headVel * dt;
            h.shankPos += h.shankVel * dt;
        }
        else
        {
            h.velocity -= hammerForce / key.hammerMass * dt;
            h.position += h.velocity * dt;
        }
    }

    // Force on the bridge from this key's strings, and the small sideways push the hammer
    // gives each string (it never strikes perfectly square)
    for (int i = 0; i < ks.numStrings; ++i)
    {
        const auto& s = strings[(size_t) (ks.firstString + i)];
        double f = 0;
        forModes (s, s.count, [&] (size_t j) { f += m.fk[j] * m.zi[j]; });
        forceSum += f;
        forceLeft += f * s.panLeft;
        forceRight += f * s.panRight;

        if (s.horizontal >= 0 && h.forceSum[(size_t) i] > 0)
        {
            const auto& hs = strings[(size_t) s.horizontal];
            const float acc = (float) (hs.hammerShare * h.forceSum[(size_t) i] / os / hs.mass);
            forModes (hs, hs.count, [&] (size_t j)
            {
                m.zr[j] += m.gr[j] * m.phi[j] * acc;
                m.zi[j] += m.gi[j] * m.phi[j] * acc;
            });
        }
        h.forceSum[(size_t) i] = 0;
    }

    // The hammer is done once it has been clear of the strings and falling for 2 ms
    ++h.samplesSinceStrike;
    h.samplesApart = (! touching && h.velocity < 0) ? h.samplesApart + 1 : 0;
    if (h.samplesApart > (int) (0.002 * sampleRate) || h.samplesSinceStrike > (int) (0.12 * sampleRate))
    {
        h.active = false;
        checkActivity();
    }
}

// While a loud note rings, the string is stretched and its tension raised. Every partial
// rises in pitch by half the relative tension rise (f ~ sqrt(T)), so a hard strike starts a
// little sharp and glides down. The stretched string also sheds its high partials faster.
void ModalPiano::applyTension()
{
    const double glide = settings.tensionGlide;
    const double loss = settings.tensionLoss;
    for (int index : longList)
    {
        const auto& head = strings[(size_t) index];
        auto& ks = keyState[(size_t) head.key];
        const double rise = settings.longitudinal ? head.tensionRise : 0.0;
        if (rise < 1e-6 && ! ks.tensionApplied)
            continue;
        ks.tensionApplied = rise >= 1e-6;
        const float shift = (float) (0.5 * rise * glide * Ts);
        const float decay = (float) (rise * loss * Ts);
        forEachString (ks, [&] (String& s)
        {
            const auto& prSrc = s.damped ? m.prDamp : m.prFree;
            const auto& piSrc = s.damped ? m.piDamp : m.piFree;
            forModes (s, s.count, [&] (size_t j)
            {
                const float w = m.omega[j];
                const float y = w * shift;
                const float a = 1.0f - w * decay - 0.5f * y * y;
                m.pr[j] = a * prSrc[j] - y * piSrc[j];
                m.pi[j] = a * piSrc[j] + y * prSrc[j];
            });
        });
    }
}

double ModalPiano::runLongitudinal (const String& s)
{
    // u_k = k y_k, so the string's squared slope is a sum of products u_k u_n
    const int n = s.longPartials;
    float* __restrict u = scratchU.data();
    float squares = 0;
    for (int k = 0; k < n; ++k)
    {
        u[k] = (float) (k + 1) * m.zi[(size_t) modeIndex (s, k)];
        squares += u[k] * u[k];
    }

    // Slowly varying part: the average tension rise, the same all along the string
    double force = s.tensionCoef * squares;
    s.tensionRise = (float) (force / s.restTension);

    // Once the stretch is tiny (below about -90 dB of the string's own force), the phantom
    // tones are inaudible: let the lengthwise modes rest and skip the pair sums
    if (s.tensionRise < 1e-7f)
    {
        for (int i = 0; i < s.longCount; ++i)
            lm.zr[(size_t) (s.longFirst + i)] = lm.zi[(size_t) (s.longFirst + i)] = 0;
        return force * settings.longitudinalCoupling * keys[(size_t) s.key].numStrings * s.outScale;
    }

    double resonant = 0;
    // Resonant part: each longitudinal mode m is driven by partial pairs whose numbers
    // differ by m (and by pairs that add up to m)
    for (int i = 0; i < s.longCount; ++i)
    {
        const size_t j = (size_t) (s.longFirst + i);
        const int order = lm.order[j];
        float pairs = 2.0f * lagProduct (u, n, order);
        for (int k = 1; k < order && k < n && order - k <= n; ++k)
            pairs += u[k - 1] * u[order - k - 1];
        const float drive = lm.drive[j] * pairs;

        const float r = lm.zr[j], im = lm.zi[j];
        const float ni = lm.pr[j] * im + lm.pi[j] * r + lm.gi[j] * drive;
        lm.zr[j] = lm.pr[j] * r - lm.pi[j] * im + lm.gr[j] * drive;
        lm.zi[j] = ni;
        resonant += lm.out[j] * (ni - drive * lm.invW2[j]);
    }
    // The lengthwise resonances reach the board less well than the slow stretch: they are
    // sharp peaks in the treble region and, coupled fully, ring as a whistle
    force += settings.longitudinalResonantShare * resonant;
    return force * settings.longitudinalCoupling * keys[(size_t) s.key].numStrings * s.outScale;
}

void ModalPiano::process (float* left, float* right, int numSamples)
{
    // Dampers of released keys land once the key has risen far enough
    for (size_t k = 0; k < keyState.size(); ++k)
    {
        auto& ks = keyState[k];
        if (ks.damperPending <= 0)
            continue;
        ks.damperPending -= numSamples;
        if (ks.damperPending <= 0)
        {
            ks.damperPending = 0;
            if (keys[k].hasDamper && ! sustain && ! ks.sostenuto)
            {
                queueAction (2, 0, 0.3, (int) k);
                // How much the strings are still moving decides how loud the release is
                double e = 0;
                for (int i = 0; i < keyState[k].numStrings; ++i)
                {
                    const auto& st = strings[(size_t) (keyState[k].firstString + i)];
                    for (int j = 0; j < st.sympatheticCount; ++j)
                    {
                        const float zr = m.zr[(size_t) (st.first + j)], zi = m.zi[(size_t) (st.first + j)];
                        e += zr * zr + zi * zi;
                    }
                }
                const double motion = std::min (2.0, std::sqrt (e) / 1.0e-3);
                if (motion > 0.01)
                    queueAction (3, 0, 0.35 * motion, (int) k);
            }
            updateDamping ((int) k);
            checkActivity();
        }
    }
    sampleClock += numSamples;

    float* __restrict zr = m.zr.data();
    float* __restrict zi = m.zi.data();
    const float* __restrict pr = m.pr.data();
    const float* __restrict pi = m.pi.data();
    const float* __restrict bgr = m.bgr.data();
    const float* __restrict bgi = m.bgi.data();
    const float* __restrict fk = m.fk.data();

    const size_t boardModes = board.omega.size();

    for (int n = 0; n < numSamples; ++n)
    {
        const float ah = (float) horizontalAcceleration;
        zoneForce.fill (0);
        double forceSum = 0, forceLeft = 0, forceRight = 0;  // vertical
        double sideSum = 0, sideLeft = 0, sideRight = 0;     // horizontal

        // Strings without a hammer on them: one step per sample
        for (const auto& r : runs)
        {
            const int j = r.first;
            const float f = stepModes (zr + j, zi + j, pr + j, pi + j, bgr + j, bgi + j, fk + j, r.count,
                                       r.side ? ah : (float) zoneAcceleration[(size_t) r.zone]);
            if (r.side)
            {
                sideSum += f;
                sideLeft += f * r.panLeft;
                sideRight += f * r.panRight;
            }
            else
            {
                forceSum += f;
                forceLeft += f * r.panLeft;
                forceRight += f * r.panRight;
                zoneForce[(size_t) r.zone] += f;
            }
        }

        // Strings under a hammer: oversampled
        for (size_t i = 0; i < hammerList.size(); ++i)
        {
            const int z = keyZone[(size_t) hammerList[i]];
            const double before = forceSum;
            runHammerKey (hammerList[i], zoneAcceleration[(size_t) z], forceSum, forceLeft, forceRight);
            zoneForce[(size_t) z] += forceSum - before;
        }

        if (++tensionCounter >= 16)
        {
            tensionCounter = 0;
            applyTension();
        }

        // Lengthwise force, through the bridge
        for (int index : longList)
        {
            const auto& s = strings[(size_t) index];
            const double f = runLongitudinal (s);
            forceSum += f;
            zoneForce[(size_t) s.zone] += f;

            // Large-motion nonlinearity (Kirchhoff-Carrier): the stretch raises the tension
            // instant by instant, and every partial feels it as a stiffer spring,
            // a_n = -w_n^2 (dT/T) y_n. Because dT itself swings at the sum and difference
            // of the partials' frequencies, this mixes them: loud notes gain a faint halo of
            // intermodulation and their pitch rises and settles as the stretch decays.
            // Updated every other sample (its effect is a slow modulation), and only while the
            // note is loud enough for it to matter
            const double rise = 4.0 * settings.nonlinear * s.tensionRise;
            if ((sampleClock + n) % 4 == 0 && rise > 4e-6)
            {
                const auto& ks = keyState[(size_t) s.key];
                const float r = (float) rise;
                const float* __restrict om = m.omega.data();
                const float* __restrict gr = m.gr.data();
                const float* __restrict gi = m.gi.data();
                auto apply = [&] (int first, int count)
                {
                    for (int j = first; j < first + count; ++j)
                    {
                        const float acc = -om[j] * om[j] * r * zi[j];
                        zr[j] += gr[j] * acc;
                        zi[j] += gi[j] * acc;
                    }
                };
                forEachStringConst (ks, [&] (const String& st)
                {
                    apply (st.first, st.sympatheticCount);
                    if (st.activity == Activity::Full)
                        apply (st.highFirst, st.count - st.sympatheticCount);
                });
            }
            forceLeft += f * s.panLeft;
            forceRight += f * s.panRight;
        }

        // The far ends of the strings, on the iron frame
        double frameLeft = 0, frameRight = 0;
        if (settings.farEnd > 0)
            runFrame (frameLeft, frameRight);

        // Forces that reach the board through the frame rather than one spot on the bridge:
        // the action's knocks, and the air under the board pushing back on it
        {
            const double before = forceSum;
            runAction (forceSum, forceLeft, forceRight);
            forceSum -= airLoad * airPressure;
            const double spread = (forceSum - before) / kZones;
            for (auto& z : zoneForce)
                z += spread;
        }

        // Soundboard and frame modes. Each mode is pushed by every zone's force times its
        // shape there, and moves every zone by the same shape.
        double radiatedLeft = 0, radiatedRight = 0, volumeVelocity = 0;
        zoneBoardVelocity.fill (0);
        {
            // Laid out zone by zone (phi[z * modes + i]) so each loop runs over the modes
            const int nb = (int) boardModes;
            float* __restrict F = board.scratchF.data();
            float* __restrict V = board.scratchV.data();
            const float* __restrict phi = board.phi.data();
            for (int i = 0; i < nb; ++i)
                F[i] = 0;
            for (int z = 0; z < kZones; ++z)
            {
                const float fz = (float) zoneForce[(size_t) z];
                const float* __restrict ph = phi + (size_t) z * (size_t) nb;
                for (int i = 0; i < nb; ++i)
                    F[i] += ph[i] * fz;
            }
            float* __restrict wr = board.wr.data();
            float* __restrict wi = board.wi.data();
            const float* __restrict qr = board.qr.data();
            const float* __restrict qi = board.qi.data();
            const float* __restrict hr = board.hr.data();
            const float* __restrict hi = board.hi.data();
            const float* __restrict om = board.omega.data();
            const float* __restrict sg = board.sigma.data();
            const float* __restrict rl = board.radLeft.data();
            const float* __restrict rrt = board.radRight.data();
            const float* __restrict bv = board.volume.data();
            float l = 0, r = 0, vol = 0;
            for (int i = 0; i < nb; ++i)
            {
                const float x = wr[i], y = wi[i];
                const float nr = qr[i] * x - qi[i] * y + hr[i] * F[i];
                const float ni = qr[i] * y + qi[i] * x + hi[i] * F[i];
                wr[i] = nr;
                wi[i] = ni;
                const float vel = om[i] * nr - sg[i] * ni;
                V[i] = vel;
                l += rl[i] * vel;
                r += rrt[i] * vel;
                vol += bv[i] * vel;
            }
            for (int z = 0; z < kZones; ++z)
            {
                const float* __restrict ph = phi + (size_t) z * (size_t) nb;
                float acc = 0;
                for (int i = 0; i < nb; ++i)
                    acc += ph[i] * V[i];
                zoneBoardVelocity[(size_t) z] = acc;
            }
            radiatedLeft = l;
            radiatedRight = r;
            volumeVelocity = vol;
        }

        // Bridge motion, fed back to every string on the next sample. The dense, resistive
        // part of the board answers half locally (the zone's own force) and half as a whole.
        // Vertical and sideways motion are linked, which is how the slow sideways vibration
        // becomes audible.
        double weightedForce = 0;
        for (int z = 0; z < kZones; ++z)
            weightedForce += zoneCoupling[(size_t) z] * zoneForce[(size_t) z];
        for (int z = 0; z < kZones; ++z)
        {
            // The resistive part goes through the bridge's own mass (a one-pole roll-off)
            const double c = zoneCoupling[(size_t) z];
            const double res = zoneResistance[(size_t) z] * bridgeResistance * c * (0.5 * c * zoneForce[(size_t) z] + 0.5 * weightedForce);
            auto& ms = bridgeMassState[(size_t) z];
            ms = res + bridgeMassCoeff[(size_t) z] * (ms - res);
            const double v = ms + crossResistance * sideSum + zoneBoardVelocity[(size_t) z];
            zoneAcceleration[(size_t) z] = (v - lastZoneVelocity[(size_t) z]) * sampleRate;
            lastZoneVelocity[(size_t) z] = v;
        }
        const double sideVelocity = horizontalResistance * sideSum + crossResistance * forceSum;
        horizontalAcceleration = (sideVelocity - lastHorizontalVelocity) * sampleRate;
        lastHorizontalVelocity = sideVelocity;

        // Sound: soundboard modes plus the direct (dense high-frequency) bridge path
        const double directLeft = forceLeft * bridgeResistance + sideLeft * crossResistance;
        const double directRight = forceRight * bridgeResistance + sideRight * crossResistance;
        float dl = (float) directLeft, dr = (float) directRight;
        for (auto& d : diffLeft) dl = d.process (dl);
        for (auto& d : diffRight) dr = d.process (dr);
        for (auto& b : highBoardLeft) dl = b.process (dl);
        for (auto& b : highBoardRight) dr = b.process (dr);
        resistiveLpLeft = dl + lpCoeff * (resistiveLpLeft - dl);
        resistiveLpRight = dr + lpCoeff * (resistiveLpRight - dr);
        // Air under the board: driven by what the board pushes out, heard in both ears,
        // and pressing back on the board through the bridge on the next sample
        double air = 0;
        if (settings.air > 0)
        {
            const float drive = (float) volumeVelocity;
            for (auto& a : airModes)
            {
                const float y = a.b0 * (drive - a.x2) - a.a1 * a.y1 - a.a2 * a.y2;
                a.x2 = a.x1; a.x1 = drive; a.y2 = a.y1; a.y1 = y;
                air += a.gain * y;
            }
            air *= settings.air;
        }
        airPressure = air;

        double outL = radiatedLeft + resistiveLpLeft + air + frameLeft;
        double outR = radiatedRight + resistiveLpRight + 0.9 * air + frameRight;

        // Directivity: the high part beams out past the open side
        beamLpLeft = outL + beamCoeff * (beamLpLeft - outL);
        beamLpRight = outR + beamCoeff * (beamLpRight - outR);
        outL -= 0.2 * (outL - beamLpLeft);
        outR += 0.25 * (outR - beamLpRight);

        // Lid reflection
        if (settings.lid > 0)
        {
            const int size = (int) lidBufLeft.size();
            lidBufLeft[(size_t) lidPos] = (float) outL;
            lidBufRight[(size_t) lidPos] = (float) outR;
            const float dl = lidBufLeft[(size_t) ((lidPos - lidDelayLeft + size) % size)];
            const float dr = lidBufRight[(size_t) ((lidPos - lidDelayRight + size) % size)];
            lidPos = (lidPos + 1) % size;
            lidHpLeft = lidHpCoeff * (lidHpLeft + dl - lidHpInLeft); lidHpInLeft = dl;
            lidHpRight = lidHpCoeff * (lidHpRight + dr - lidHpInRight); lidHpInRight = dr;
            lidLpLeft = lidHpLeft + lidLpCoeff * (lidLpLeft - lidHpLeft);
            lidLpRight = lidHpRight + lidLpCoeff * (lidLpRight - lidHpRight);
            // The open lid sends the board's sound towards the listener: the reflection adds
            // to the direct sound, more so on the open (right) side
            outL += 0.2 * settings.lid * lidLpRight;
            outR += 0.28 * settings.lid * lidLpLeft;
            // The lid panel's own ringing, heard in both ears
            const float drive = (float) (dl + dr);
            float panel = 0;
            for (size_t i = 0; i < lidPanel.size(); ++i)
                panel += lidPanelGain[i] * lidPanel[i].process (drive);
            outL += settings.lid * panel;
            outR += 0.8 * settings.lid * panel;
        }

        // Room: early reflections, heard from the player's seat
        if (room > 0.0f && ! reflections.empty())
        {
            const int size = (int) roomBufLeft.size();
            roomBufLeft[(size_t) roomPos] = (float) outL;
            roomBufRight[(size_t) roomPos] = (float) outR;
            float rl = 0, rr = 0;
            for (auto& r : reflections)
            {
                const int at = (roomPos - r.delay + size) % size;
                const float x = r.fromRight ? roomBufRight[(size_t) at] : roomBufLeft[(size_t) at];
                r.state = x + r.lpCoeff * (r.state - x);
                (r.toRight ? rr : rl) += r.gain * r.state;
            }
            roomPos = (roomPos + 1) % size;
            outL += room * rl;
            outR += room * rr;
        }

        // DC blocker
        dcLeft = outL - dcInLeft + dcCoeff * dcLeft;
        dcInLeft = outL;
        dcRight = outR - dcInRight + dcCoeff * dcRight;
        dcInRight = outR;

        left[n] = (float) (dcLeft * outputGain);
        right[n] = (float) (dcRight * outputGain);

        if (++activityCounter >= 256)
        {
            activityCounter = 0;
            checkActivity();
        }
    }
}

} // namespace openpiano::piano
