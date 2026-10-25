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
        const double b1 = 0.15 + 0.35 * octavesAboveA0 / 7.25;
        // Air and internal friction take the upper partials much faster (Bank 2000; Chaigne & Askenfelt 1994)
        const double b2 = wound ? 4.0e-7 : 2.0e-7;
        // Above about 2 kHz the losses climb steeply, so treble overtones die within a second
        // while the fundamental sings on
        // (the top fundamentals are spared: their losses are mostly at the bridge)
        const double overtone = std::clamp (f / f1 - 1.0, 0.0, 1.0);
        return b1 + b2 * f * f + 0.6 * overtone * std::pow (f / 2000.0, 3);
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

    // One sample of the soundboard modes driven by force F. Returns their summed velocity
    // and the velocity weighted towards each ear.
    void stepBoard (float* __restrict wr, float* __restrict wi, const float* __restrict qr,
                    const float* __restrict qi, const float* __restrict hr, const float* __restrict hi,
                    const float* __restrict omega, const float* __restrict sigma,
                    const float* __restrict radL, const float* __restrict radR, int n, float F,
                    double& velocity, double& left, double& right)
    {
        float v = 0, l = 0, r = 0;
        for (int i = 0; i < n; ++i)
        {
            const float a = wr[i], b = wi[i];
            const float nr = qr[i] * a - qi[i] * b + hr[i] * F;
            const float ni = qr[i] * b + qi[i] * a + hi[i] * F;
            wr[i] = nr;
            wi[i] = ni;
            const float vel = omega[i] * nr - sigma[i] * ni;
            v += vel;
            l += radL[i] * vel;
            r += radR[i] * vel;
        }
        velocity = v;
        left = l;
        right = r;
    }
}

void ModalPiano::Modes::clear()
{
    for (auto* v : { &zr, &zi, &pr, &pi, &prFree, &piFree, &prDamp, &piDamp, &gr, &gi, &bgr, &bgi, &fk, &beta,
                     &phi, &prOs, &piOs, &gOsr, &gOsi, &sigmaInternal, &omega })
        v->clear();
}

ModalPiano::ModalPiano() = default;

void ModalPiano::prepare (double newSampleRate, const Settings& newSettings)
{
    sampleRate = newSampleRate;
    Ts = 1.0 / sampleRate;
    settings = newSettings;
    keys = buildConcertGrandScale (settings.unisonSpreadCents);

    buildModes();
    buildSoundboard();
    buildAction();

    // Output filters: gentle top-end roll-off on the direct bridge path, and a DC blocker
    lpCoeff = std::exp (-2.0 * M_PI * 7000.0 / sampleRate);
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
    const double strike = strikeRatio * (1.0 + 0.02 * rng.symmetric());
    const double B = spec.inharmonicity * (1.0 + 0.03 * rng.symmetric());
    const double lossScale = 1.0 + 0.1 * rng.symmetric();
    const double fFlex = spec.fundamental / std::sqrt (1.0 + B);
    const double L = spec.length;

    String s;
    s.key = keyIndex;
    s.kind = kind;
    s.first = (int) m.size();
    s.mass = (float) (spec.linearDensity * L / 2.0);

    // Pan with the keyboard, as heard from the player's seat
    const double position = (key.midiNote - kLowestKey) / double (kHighestKey - kLowestKey);
    const double angle = M_PI / 4 + (position - 0.5) * 0.9 + 0.03 * rng.symmetric();
    s.panLeft = (float) std::cos (angle);
    s.panRight = (float) std::sin (angle);

    int sympathetic = 0;
    for (int n = 1; n <= maxModes; ++n)
    {
        const double f = n * fFlex * std::sqrt (1.0 + B * n * n);
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
        m.phi.push_back ((float) (std::sin (n * M_PI * strike) * width));
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
    const int sympatheticCap = kind == Kind::Duplex || ! key.hasDamper ? s.count : kind == Kind::Horizontal ? 6 : 16;
    s.sympatheticCount = ! key.hasDamper || kind == Kind::Duplex ? s.count
                                                                 : std::clamp (sympathetic, 1, std::min (s.count, sympatheticCap));

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
    const double fTop = std::min (0.4 * sampleRate, 12000.0);

    int partials = 0;
    for (int j = 0; j < s.count && j < 32; ++j)
        if (m.omega[(size_t) modeIndex (s, j)] < 2.0 * M_PI * 4000.0)
            partials = j + 1;

    s.longFirst = (int) lm.zr.size();
    s.longCount = 0;
    s.longPartials = partials;
    s.tensionCoef = (float) (EA / 4.0 * std::pow (M_PI / L, 2));
    s.restTension = (float) spec.tension;

    for (int order = 1; order <= 4 && order * f1 < fTop && partials > 1; ++order)
    {
        const double w = 2.0 * M_PI * order * f1;
        // Free longitudinal ringing dies fast (end losses at the bridge pin and agraffe); what
        // remains audible is mostly the forced, phantom part that follows the transverse motion
        const double sigma = 40.0;
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
    scratchU.resize (std::max (scratchU.size(), (size_t) 32));
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
                strings[(size_t) index].hammerShare = (float) (0.03 + 0.04 * rng.next());
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
                     &m.bgi, &m.fk, &m.beta, &m.phi, &m.prOs, &m.piOs, &m.gOsr, &m.gOsi, &m.sigmaInternal, &m.omega })
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

void ModalPiano::buildSoundboard()
{
    // A grand soundboard has a roughly constant modal density (about one mode every 30 Hz).
    // Each mode is a damped resonator driven by the total bridge force. Mode masses are set
    // so that, averaged over frequency, the modes add about 1.5e-4 m/s/N to the bridge's
    // admittance (with the resistive part, about 3e-4 in all), in line with measured bridge
    // admittances (Wogram 1980; Conklin 1996; Giordano 1998). The low modes are sparse and
    // stand out as resonances, up to 1.5e-3 m/s/N. The crown and the strings' downbearing
    // stiffen the board (below).
    Rng rng (1234);
    struct Mode { double f, loss, peak, radiation; };
    std::vector<Mode> modes;

    const double lossFactor = 0.02;
    const double meanBoardAdmittance = 0.9e-4;
    const double meanSpacing = 2.0 * M_PI * 33.0; // rad/s between modes
    // Crown and downbearing: the arched board, pressed down by the strings, is stiffer than
    // a flat panel. Stiffness adds to the mode's spring, so each frequency becomes
    // sqrt(f_flat^2 + f_c^2). Low modes rise a lot (the lowest from about 32 Hz to about 58 Hz,
    // as measured on crowned grand boards), upper ones hardly move, and the low modes spread
    // further apart.
    const double fc = 48.0 * settings.crown;
    for (double f = 32.0; f < std::min (7000.0, 0.45 * sampleRate);)
    {
        const double flat = f * (1.0 + 0.12 * (rng.next() - 0.5));
        const double fm = std::sqrt (flat * flat + fc * fc);
        f += 33.0 * (0.6 + 0.8 * rng.next());
        if (fm > 0.45 * sampleRate)
            break;
        // Same mass, stiffer spring: the mode answers a push less readily, by flat / fm
        const double sigmaFlat = M_PI * flat * lossFactor;
        const double peak = std::min (0.9e-3, meanBoardAdmittance * meanSpacing / (M_PI * sigmaFlat)) * flat / fm;
        modes.push_back ({ fm, lossFactor, peak, 1.0 });
    }

    // The rim and frame: a few heavy, lightly damped low resonances that the bridge
    // barely moves, but which keep humming under the sound
    for (double f : { 68.0, 104.0, 151.0, 207.0, 286.0, 410.0 })
        modes.push_back ({ f * (1.0 + 0.05 * rng.symmetric()), 0.008, 0.9e-4, 0.5 });

    const size_t n = modes.size();
    for (auto* v : { &board.wr, &board.wi, &board.qr, &board.qi, &board.hr, &board.hi, &board.omega, &board.sigma,
                     &board.radLeft, &board.radRight })
        v->assign (n, 0.0f);

    for (size_t i = 0; i < n; ++i)
    {
        const double f = modes[i].f;
        const double w = 2.0 * M_PI * f;
        const double sigma = M_PI * f * modes[i].loss;
        const double modalMass = 1.0 / (2.0 * sigma * modes[i].peak);

        const cd pole = std::exp (cd (-sigma, w) * Ts);
        const cd gain = pole * (Ts / (w * modalMass));
        board.omega[i] = (float) w;
        board.sigma[i] = (float) sigma;
        board.qr[i] = (float) pole.real();
        board.qi[i] = (float) pole.imag();
        board.hr[i] = (float) gain.real();
        board.hi[i] = (float) gain.imag();

        // Direction. Below a few hundred Hz the board is small next to the wavelength and
        // pushes air like one piston: both ears hear the same thing. Higher up, each mode
        // splits the board into patches moving in opposite directions, so its sound comes out
        // in lobes and differs from ear to ear. The bass end of the board is nearer the left
        // ear, the treble end nearer the right (player's view).
        const double spread = std::clamp ((f - 150.0) / 1200.0, 0.0, 1.0);
        const double common = (0.7 + 0.3 * rng.next()) * (1.0 - spread);
        const double sideL = (0.4 + 0.6 * rng.next()) * (rng.next() < 0.5 ? -1 : 1) * spread;
        const double sideR = (0.4 + 0.6 * rng.next()) * (rng.next() < 0.5 ? -1 : 1) * spread;
        const double tilt = 0.15 * std::clamp ((f - 100.0) / 2000.0, -1.0, 1.0);
        const double l = (common + sideL) * (1.0 - tilt);
        const double r = (common + sideR) * (1.0 + tilt);
        board.radLeft[i] = (float) (l * modes[i].radiation);
        board.radRight[i] = (float) (r * modes[i].radiation);
    }

    // Overall level: a mezzo-forte middle C peaks around -20 dBFS, fortissimo around -10
    outputGain = 50.0;
    buildBody();
}

void ModalPiano::buildBody()
{
    // Air under a grand's board is boxed in by the rim and the floor and leaks out at the
    // gap underneath: a few broad resonances between about 90 and 350 Hz that fill out the bass.
    airModes.clear();
    const double airF[] = { 92.0, 165.0, 238.0, 330.0 };
    const double airQ[] = { 4.0, 5.0, 5.0, 6.0 };
    const double airG[] = { 0.55, 0.4, 0.3, 0.2 };
    for (int i = 0; i < 4; ++i)
    {
        const double w = 2.0 * M_PI * airF[i] / sampleRate;
        const double alpha = std::sin (w) / (2.0 * airQ[i]);
        const double a0 = 1.0 + alpha;
        airModes.push_back ({ (float) (alpha / a0), (float) (-2.0 * std::cos (w) / a0), (float) ((1.0 - alpha) / a0),
                              0, 0, 0, 0, (float) airG[i] });
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
}

void ModalPiano::buildAction()
{
    // Wooden action parts ring briefly at a few hundred Hz to a few kHz (Askenfelt 1993,
    // "Observations on the transient components of the piano tone"). Three kinds of hit:
    //   0 key bottom: the key lands on the keybed felt, a dull low thump
    //   1 hammer knock: shank and hammer flex at the strike, a brighter click
    //   2 damper return: the felt lands back on the strings, soft and muffled
    const double freqs[8] = { 95, 160, 270, 430, 690, 1100, 1750, 2800 };
    const double q[8] = { 6, 7, 8, 9, 10, 10, 11, 12 };
    const double shapes[3][8] = {
        { 1.0, 0.9, 0.7, 0.45, 0.25, 0.12, 0.05, 0.02 },
        { 0.15, 0.25, 0.4, 0.6, 0.8, 0.9, 0.7, 0.5 },
        { 0.4, 0.5, 0.45, 0.3, 0.15, 0.06, 0.02, 0.01 },
    };
    for (int i = 0; i < 8; ++i)
    {
        const double w = 2.0 * M_PI * freqs[i];
        const cd pole = std::exp (cd (-w / (2.0 * q[i]), w) * Ts);
        action.pr[(size_t) i] = (float) pole.real();
        action.pi[(size_t) i] = (float) pole.imag();
        for (int t = 0; t < 3; ++t)
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
        for (size_t j = 0; j < 8; ++j)
            action.zr[j] += h.amount * action.weight[(size_t) h.type][j];
        // The bank radiates from roughly where the keys were hit
        action.panLeft += 0.5f * (h.panLeft - action.panLeft);
        action.panRight += 0.5f * (h.panRight - action.panRight);
        h = actionHits.back();
        actionHits.pop_back();
    }
    float f = 0;
    for (size_t j = 0; j < 8; ++j)
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
    lastHorizontalVelocity = horizontalAcceleration = 0;
    resistiveLpLeft = resistiveLpRight = 0;
    buildBody();
    dcLeft = dcRight = dcInLeft = dcInRight = 0;
    sustain = false;

    for (auto& k : keyState)
    {
        k.down = false;
        k.hammer = {};
    }
    for (auto& s : strings)
    {
        s.activity = Activity::Off;
        s.damped = true;
    }
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
    const bool damped = key.hasDamper && ! ks.down && ! sustain;

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
    // Releasing a key: the damper drops back onto the strings (if it was off them)
    if (! down && keyState[(size_t) k].down && keys[(size_t) k].hasDamper && ! sustain)
        queueAction (2, (int) (0.01 * sampleRate), 0.3, k);
    keyState[(size_t) k].down = down;
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

void ModalPiano::strike (int midiNote, double hammerVelocity)
{
    const int k = midiNote - kLowestKey;
    if (k < 0 || k >= (int) keys.size())
        return;

    auto& ks = keyState[(size_t) k];
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
        queueAction (1, 0, 0.06 * std::pow (v, 1.5) * (0.4 + 0.6 * treble), k);
    }

    // Soft notes are made a little louder without changing their tone: the key's strings push
    // harder on the bridge and feel its motion less by the same factor, so their decay and
    // coupling stay as they were. Fortissimo (5 m/s) is left unchanged.
    {
        const float lift = (float) std::pow (5.0 / std::clamp (hammerVelocity, 0.3, 5.0), settings.softLift);
        forEachString (ks, [&] (String& s)
        {
            const float ratio = lift / s.outScale;
            forModes (s, s.count, [&] (size_t j) { m.fk[j] *= ratio; m.bgr[j] /= ratio; m.bgi[j] /= ratio; });
            s.outScale = lift;
        });
    }

    ks.hammer.active = true;
    ks.hammer.position = nearest; // touching the first string
    ks.hammer.velocity = std::max (0.05, hammerVelocity);
    ks.hammer.headPos = ks.hammer.shankPos = nearest;
    ks.hammer.headVel = ks.hammer.shankVel = ks.hammer.velocity;
    ks.hammer.samplesSinceStrike = 0;
    ks.hammer.samplesApart = 0;
    for (int i = 0; i < ks.numStrings; ++i)
    {
        const auto& s = strings[(size_t) (ks.firstString + i)];
        double u = 0;
        forModes (s, s.count, [&] (size_t j) { u += m.phi[j] * m.zi[j]; });
        ks.hammer.lastCompression[(size_t) i] = nearest - u - ks.hammer.gap[(size_t) i];
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
        fullList.push_back (i);

    // Merge the computed mode ranges into runs where they sit back to back, at most 12
    // strings (about 4 keys) each, so their pan barely differs
    auto& segments = segmentScratch;
    segments.clear();
    for (int index : fullList)
    {
        const auto& s = strings[(size_t) index];
        const bool side = s.kind == Kind::Horizontal;
        segments.push_back ({ s.first, s.sympatheticCount, s.panLeft, s.panRight, side });
        if (s.activity == Activity::Full && s.count > s.sympatheticCount)
            segments.push_back ({ s.highFirst, s.count - s.sympatheticCount, s.panLeft, s.panRight, side });
    }
    std::sort (segments.begin(), segments.end(), [] (const Segment& a, const Segment& b) { return a.first < b.first; });

    runs.clear();
    int runSegments = 0;
    for (const auto& g : segments)
    {
        if (! runs.empty() && runs.back().first + runs.back().count == g.first && runs.back().side == g.side
            && runSegments < 12)
        {
            auto& r = runs.back();
            r.count += g.count;
            r.panLeft += (g.panLeft - r.panLeft) / float (runSegments + 1);
            r.panRight += (g.panRight - r.panRight) / float (runSegments + 1);
            ++runSegments;
        }
        else
        {
            runs.push_back ({ g.first, g.count, g.panLeft, g.panRight, g.side });
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
    const double q0 = key.feltStiffness * hardness;
    const double p = key.feltExponent;
    // Felt has memory (Stulov 1995): squeezed quickly it is stiff; squeezed slowly it has
    // time to relax and gives way. So a gentle touch makes a long, soft contact that only
    // excites the low partials, and a hard blow a short, stiff one that brings out the top.
    const double memory = settings.feltMemory;
    const double coreDepth = key.outerFeltDepth * 1000.0 * settings.outerFeltScale; // [mm]
    // Treble hammers are smaller and their outer felt firmer, or soft notes there would vanish
    const double outer = settings.outerFelt * std::pow (10.0, std::max (0.0, (keyIndex - 39) / 48.0) * trebleFirmness);
    const double relax = 1.0 - std::exp (-dt / std::max (1e-6, settings.feltRelaxTime));

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
            const double f = compression > 0 ? std::max (0.0, q0 * (elastic - memory * history)) : 0.0;
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
    if (h.samplesApart > (int) (0.002 * sampleRate) || h.samplesSinceStrike > (int) (0.1 * sampleRate))
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
        force += lm.out[j] * (ni - drive * lm.invW2[j]);
    }
    return force * settings.longitudinalCoupling * keys[(size_t) s.key].numStrings * s.outScale;
}

void ModalPiano::process (float* left, float* right, int numSamples)
{
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
        const float av = (float) bridgeAcceleration;
        const float ah = (float) horizontalAcceleration;
        double forceSum = 0, forceLeft = 0, forceRight = 0;  // vertical
        double sideSum = 0, sideLeft = 0, sideRight = 0;     // horizontal

        // Strings without a hammer on them: one step per sample
        for (const auto& r : runs)
        {
            const int j = r.first;
            const float f = stepModes (zr + j, zi + j, pr + j, pi + j, bgr + j, bgi + j, fk + j, r.count, r.side ? ah : av);
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
            }
        }

        // Strings under a hammer: oversampled
        for (size_t i = 0; i < hammerList.size(); ++i)
            runHammerKey (hammerList[i], bridgeAcceleration, forceSum, forceLeft, forceRight);

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
            forceLeft += f * s.panLeft;
            forceRight += f * s.panRight;
        }

        runAction (forceSum, forceLeft, forceRight);

        // The air under the board pushes back on it
        forceSum -= airLoad * airPressure;

        // Soundboard
        double boardVelocity = 0, radiatedLeft = 0, radiatedRight = 0;
        stepBoard (board.wr.data(), board.wi.data(), board.qr.data(), board.qi.data(), board.hr.data(),
                   board.hi.data(), board.omega.data(), board.sigma.data(), board.radLeft.data(),
                   board.radRight.data(), (int) boardModes, (float) forceSum, boardVelocity, radiatedLeft,
                   radiatedRight);

        // Bridge motion, fed back to every string on the next sample. Vertical and sideways
        // motion are linked, which is how the slow sideways vibration becomes audible.
        const double bridgeVelocity = bridgeResistance * forceSum + crossResistance * sideSum + boardVelocity;
        const double sideVelocity = horizontalResistance * sideSum + crossResistance * forceSum;
        bridgeAcceleration = (bridgeVelocity - lastBridgeVelocity) * sampleRate;
        horizontalAcceleration = (sideVelocity - lastHorizontalVelocity) * sampleRate;
        lastBridgeVelocity = bridgeVelocity;
        lastHorizontalVelocity = sideVelocity;

        // Sound: soundboard modes plus the direct (dense high-frequency) bridge path
        const double directLeft = forceLeft * bridgeResistance + sideLeft * crossResistance;
        const double directRight = forceRight * bridgeResistance + sideRight * crossResistance;
        resistiveLpLeft = directLeft + lpCoeff * (resistiveLpLeft - directLeft);
        resistiveLpRight = directRight + lpCoeff * (resistiveLpRight - directRight);
        // Air under the board: driven by what the board pushes out, heard in both ears,
        // and pressing back on the board through the bridge on the next sample
        double air = 0;
        if (settings.air > 0)
        {
            const float drive = (float) (0.5 * (radiatedLeft + radiatedRight));
            for (auto& a : airModes)
            {
                const float y = a.b0 * (drive - a.x2) - a.a1 * a.y1 - a.a2 * a.y2;
                a.x2 = a.x1; a.x1 = drive; a.y2 = a.y1; a.y1 = y;
                air += a.gain * y;
            }
            air *= settings.air;
        }
        airPressure = air;

        double outL = radiatedLeft + resistiveLpLeft + air;
        double outR = radiatedRight + resistiveLpRight + 0.9 * air;

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
            outL += 0.35 * settings.lid * lidLpRight;
            outR += 0.5 * settings.lid * lidLpLeft;
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
