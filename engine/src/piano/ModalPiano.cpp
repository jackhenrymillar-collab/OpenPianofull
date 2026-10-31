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
#include "FrameModes.inc"
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
    // String losses from first principles (Valette & Cuesta 1993; Woodhouse 2004): the air's
    // viscous drag on the wire, from its diameter and the air's density and viscosity, plus the
    // steel's internal loss: a viscoelastic loss on the bending stiffness (strongest where the
    // partials are stiff, i.e. high partials and short strings) and a small frequency-independent
    // loss (dislocation friction, and the winding's rubbing on a wound string). The two material
    // constants are the same for every key.
    constexpr double kAirDensity = 1.2, kAirViscosity = 1.5e-5;   // kg/m^3, m^2/s
    double internalDecay (double f, const StringSpec& spec)
    {
        // Material constants, chosen once for the whole keyboard on the late decay rates of the
        // octave bands (A0-B1 2.7, C6-B6 6, C7-C8 12 dB/s on the reference grand): the
        // viscoelastic loss factor of the steel and a constant loss per cycle for plain and for
        // copper-wound strings (the winding rubs on the core)
        constexpr double dVe = 2.2e-3, d0P = 1.0e-4, d0W = 3.0e-4;
        const double w = 2.0 * M_PI * f;
        const double d = spec.outerDiameter;
        const double rhoEff = spec.linearDensity / (M_PI * d * d / 4.0);              // mass per volume the air sees
        const double M = (2.0 / d) * std::sqrt (2.0 * kAirViscosity / w);             // boundary layer over the radius
        const double qAir = (kAirDensity / rhoEff) * (2.0 * std::sqrt (2.0) * M + M * M);
        const double I = M_PI * std::pow (spec.coreDiameter, 4) / 64.0;             // core bending stiffness
        const double qVe = dVe * 4.0 * M_PI * M_PI * f * f * spec.linearDensity * kSteelE * I / (spec.tension * spec.tension);
        const double q0 = spec.wound ? d0W : d0P;
        return M_PI * f * (qAir + qVe + q0);   // loss factor to decay rate [Np/s]
    }

    // The felt takes this long to settle into its full grip once the damper has landed
    constexpr int kDamperRampMs = 40;

    // Extra decay with the felt damper resting on the string. A real damper takes a note
    // 25-30 dB down within the first 0.3 s after release (Salamander's release samples); with
    // the felt settling over kDamperRampMs that needs about 70 nepers/s on the low partials.
    // 10 left the damped note ringing 10-15 dB louder than the recordings under the next notes.
    double damperDecay (double f, double strength)
    {
        return strength * 70.0 * (1.0 + f / 2000.0);
    }

    // One sample of a run of string modes driven by the bridge acceleration a.
    // Returns their total force on the bridge.
    // The same step for a string whose input this sample is not known yet: last sample's
    // acceleration is folded into the state first, and the force returned is the one without
    // this sample's input; gain is the force per unit of that input (the direct part, which the
    // bridge solve adds once the acceleration is known)
    float stepModesDeferred (float* __restrict zr, float* __restrict zi, const float* __restrict pr,
                             const float* __restrict pi, const float* __restrict bgr, const float* __restrict bgi,
                             const float* __restrict fk, int n, float a, float& gain)
    {
        float f = 0, g = 0;
        for (int j = 0; j < n; ++j)
        {
            const float r = zr[j] + bgr[j] * a, im = zi[j] + bgi[j] * a;
            const float ni = pr[j] * im + pi[j] * r;
            zr[j] = pr[j] * r - pi[j] * im;
            zi[j] = ni;
            f += fk[j] * ni;
            g += fk[j] * bgi[j];
        }
        gain = g;
        return f;
    }

    // Solves A x = b in place (A n by n, row-major; b becomes x), Gaussian elimination with
    // partial pivoting: the bridge system is small (the live zones) and well conditioned
    void solveDense (double* A, double* b, int n)
    {
        for (int k = 0; k < n; ++k)
        {
            int piv = k;
            for (int i = k + 1; i < n; ++i)
                if (std::abs (A[i * n + k]) > std::abs (A[piv * n + k]))
                    piv = i;
            if (piv != k)
            {
                for (int j = 0; j < n; ++j)
                    std::swap (A[k * n + j], A[piv * n + j]);
                std::swap (b[k], b[piv]);
            }
            const double d = A[k * n + k];
            if (d == 0.0)
                continue;
            for (int i = k + 1; i < n; ++i)
            {
                const double m = A[i * n + k] / d;
                if (m == 0.0)
                    continue;
                for (int j = k; j < n; ++j)
                    A[i * n + j] -= m * A[k * n + j];
                b[i] -= m * b[k];
            }
        }
        for (int k = n - 1; k >= 0; --k)
        {
            double v = b[k];
            for (int j = k + 1; j < n; ++j)
                v -= A[k * n + j] * b[j];
            const double d = A[k * n + k];
            b[k] = d != 0.0 ? v / d : 0.0;
        }
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
                     &phi, &pt, &wa, &prOs, &piOs, &gOsr, &gOsi, &sigmaInternal, &omega, &alt, &damperShift })
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
    // Every zone couples to the board exactly as the solved shapes say: the bridge sections are
    // the Steinway D's, so no fitted per-zone factor is left
    zoneCoupling.fill (1.0);

    // The tuner tunes the piano as built, by ear: the body pulls some notes a few cents
    // sharp or flat (most near strong body modes), and the pins are turned to cancel it.
    // Measured on this model with the pitch test (tools: openpiano-piano-test), in cents.
    tuneFactor.assign (keys.size(), 1.0);
    for (size_t k = 0; k < keys.size() && k < kTunerCents.size(); ++k)
        tuneFactor[k] = std::pow (2.0, -kTunerCents[k] / 1200.0);
    buildModes();
    for (auto& st : strings)
        st.zone = keyZone[(size_t) st.key];
    buildSoundboard();
    buildFrame();
    buildAction();

    // Output filters: gentle top-end roll-off on the direct bridge path, and a DC blocker
    lpCoeff = std::exp (-2.0 * M_PI * 3500.0 / sampleRate);
    lp2Coeff = 0.0;   // the two extra poles above 8 kHz were a fit to the C5 recordings' top band: off
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
    // percent, and stiffness goes with its fourth power: B differs by about 1% (more split the treble
    // 2nd partials into separate beating peaks the recordings don't show). That keeps
    // the upper partials of a unison from all beating in step, which would dull the tone at once.
    const double B = spec.inharmonicity * (1.0 + 0.01 * rng.symmetric());
    const double lossScale = (1.0 + 0.1 * rng.symmetric()) * (1.0 + 1.2 * ageFactor);
    const double fFlex = spec.fundamental * (keyIndex < (int) tuneFactor.size() ? tuneFactor[(size_t) keyIndex] : 1.0) / std::sqrt (1.0 + B);
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
        const double sigma = internalDecay (f, spec) * lossScale;
        const double sign = (n % 2 == 1) ? 1.0 : -1.0; // (-1)^(n+1)

        const cd pole = std::exp (cd (-sigma, w) * Ts);
        const cd poleOs = std::exp (cd (-sigma, w) * (Ts / os));
        const cd gain = pole * (Ts / w);          // impulse-invariant acceleration input
        const cd gainOs = poleOs * (Ts / os / w);
        // The bridge's acceleration each sample is the change of its velocity over that sample,
        // so it belongs to the sample's middle: the free strings (stepModesDeferred) take it
        // half a sample before the state they report. Taken a whole sample early, as the state
        // update alone would have it, the string's force against the bridge's velocity has a
        // negative real part just below each partial (0.1% of its peak at 1 kHz, 60% at 20 kHz),
        // and with the bridge's sideways mobility above 3 kHz that fed a note's own sideways
        // modes until they grew (C6 on the photo-traced body: +4 dB/s). At half a sample the
        // discrete admittance is positive real at every frequency, like the string's own.
        const cd gainBridge = std::exp (cd (-sigma, w) * (0.5 * Ts)) * (Ts / w);

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
        m.bgr.push_back ((float) (gainBridge.real() * beta)); m.bgi.push_back ((float) (gainBridge.imag() * beta));
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
        {
            // The damper resting on the string: its felt (a soft spring, kd per string) carrying
            // the damper head's mass (Md, shared by the key's strings), pressing just behind the
            // hammer line. Below the felt-and-head resonance the head's mass pulls the partials
            // down; above it the felt's spring pushes them up, towards the pitch of the length
            // behind the damper. The hard limit is that length's pitch (the felt as a node).
            const double xd = std::clamp (strike + 0.07, 0.1, 0.3);
            // (The felt is soft under the damper's own light weight: a few hundred N/m per
            // string. 3000 N/m shifted the tenor's damped partials by up to 165 cents, heard
            // as an off-pitch flick on every release.)
            const double kd = 800.0, eta = 0.6, Md = 0.012 / std::max (1, key.numStrings);
            const cd kf (kd, kd * eta);
            const cd keff = kf * (-w * w * Md) / (kf - w * w * Md);
            const double sn = std::sin (n * M_PI * xd);
            const double shift = keff.real() * sn * sn / (w * w * spec.linearDensity * L);
            m.damperShift.push_back ((float) std::clamp (shift, -0.3, 1.0 / (1.0 - xd) - 1.0));
        }

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
    s.sympTopF = m.omega[(size_t) (s.first + s.sympatheticCount - 1)] / (float) (2.0 * M_PI);

    // Voicing: a technician evens out a piano's loudness note by note (hammer needling,
    // string levelling). This curve was fitted to the recorded grand's level, note by note
    // relative to the middle of the range. Pushing harder on the bridge while feeling its
    // motion less by the same factor keeps the decay and coupling unchanged.
    {
        static constexpr double voicingDb[30] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // every 3 keys from A0
        const double pos = std::clamp ((key.midiNote - kLowestKey) / 3.0, 0.0, 28.999);
        const int i = (int) pos;
        const double db = voicingDb[i] + (voicingDb[i + 1] - voicingDb[i]) * (pos - i);
        const float g = (float) std::pow (10.0, db / 20.0);
        // Tone voicing (needling or hardening the felt): the slope of the upper partials, in dB
        // per octave of partial number, fitted note by note to the recorded grand's brightness
        static constexpr double toneDb[30] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // every 3 keys from A0
        const double tone = toneDb[i] + (toneDb[i + 1] - toneDb[i]) * (pos - i);
        // In the top octave and a half the bridge and board give way less and less to the upper
        // partials (up to 15 dB less per octave of partial number from C7 up). Without this the
        // 2nd and 3rd partials there sit 10 to 20 dB above the recorded grand's: thin and metallic.
        // The cut peaks at C7 and fades out by F#7: above there the recorded partials are no
        // stronger than ours, and cutting them leaves a bare, harp-like sine.
        const double topRoll = 0.0;   // was a 15 dB/octave cut from C7 fitted to the C5 recordings; the D's own top strings and bridge decide now
        for (int j = 0; j < s.count; ++j)
        {
            const size_t q = (size_t) (s.first + j);
            const float gj = g * (float) std::pow (10.0, (tone - topRoll) * std::log2 (j + 1.0) / 20.0);
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
                     &m.bgi, &m.fk, &m.beta, &m.phi, &m.pt, &m.wa, &m.prOs, &m.piOs, &m.gOsr, &m.gOsi, &m.sigmaInternal, &m.omega, &m.alt, &m.damperShift })
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

    for (const auto& s : strings)
    {
        // Treble dampers are smaller and lighter: the recordings' released notes sit about
        // 6 dB higher there than in the tenor, so the felt's grip tapers from C4 to C6.
        const double taper = 1.0 - 0.4 * std::clamp ((s.key - 39) / 24.0, 0.0, 1.0);
        forModes (s, s.count, [&] (size_t j)
        {
            const double w = m.omega[j] * (1.0 + m.damperShift[j]);
            const double sigma = m.sigmaInternal[j] + taper * damperDecay (w / (2.0 * M_PI), damperStrength);
            const cd pole = std::exp (cd (-sigma, w) * Ts);
            m.prDamp[j] = (float) pole.real();
            m.piDamp[j] = (float) pole.imag();
        });
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
    // The solved bridge admittance is about 25x what the recordings imply (likely the
    // bridge and string-load mass missing from the solve), so the coupling is fitted down.
    const double phiScale = 1.0;   // the solved admittance, as it is
    // The soundboard and the body around it: modes of a Steinway D's body (RT-425 geometry), computed by
    // tools/soundboard/build_board.py (BoardModes.inc). The model is a tapered, ribbed spruce
    // board with its bridges, crown and downbearing, glued into a rim that bends and twists,
    // held by back braces, standing on three legs, with the iron frame bolted on and the air
    // above and below the board, all solved together. Everything a mode does comes from that
    // model: its shape at the bridge zones (unit modal mass, so the admittance the strings
    // feel is the physical one), its losses (the materials it strains, the sound it radiates,
    // the air leaving the case) and the sound it sends to the listener. Modes that mostly move
    // the frame are in FrameModes.inc instead (see buildFrame).
    static_assert (kBoardZones == kZones, "regenerate BoardModes.inc for this many bridge zones");
    constexpr int Z = kZones;
    struct Mode { double f, loss; int n; std::array<double, Z> phi; };
    std::vector<Mode> modes;

    // The crown and downbearing are in the model's stiffness. An old board has lost some of its
    // crown, which lowers its low modes a little; the Crown setting scales what is left of that
    // stiffening (1 = as built). Age also makes the wood and glue lossier.
    const double crownLift = settings.crown * (1.0 - 0.5 * ageFactor);
    const double ageLoss = 1.0 + 0.8 * ageFactor;
    // Modes are simulated one by one up to 5 kHz; above that they overlap so much that the
    // board acts as a plain resistance (the direct bridge path below)
    const double fTop = std::min (5000.0, 0.45 * sampleRate);
    for (int i = 0; i < kBoardModeCount; ++i)
    {
        const auto& src = kBoardModes[i];
        // The crown's share of each mode's stiffness falls off with frequency (about 12 Hz equivalent)
        const double fm = std::sqrt (std::max (1.0, (double) src.f * src.f + (crownLift - 1.0) * 144.0));
        if (fm > fTop)
            continue;
        // The casters rub on the floor: the case's rocking modes, below ~60 Hz, lose far more
        // than the wood alone would give them
        const double floorLoss = 0.3 / (1.0 + std::pow (fm / 90.0, 4.0));   // the case rocks on its casters and legs: heavily damped (lighter damping let these modes pull bass notes up to 25 cents off)
        const double loss = src.radiationLoss + ageLoss * src.structLoss + src.airLoss + floorLoss;
        Mode md { fm, loss, i, {} };
        for (int z = 0; z < Z; ++z)
            md.phi[(size_t) z] = src.phi[z];
        modes.push_back (md);
    }

    std::sort (modes.begin(), modes.end(), [] (const Mode& a, const Mode& b) { return a.f < b.f; });   // (zoneModesV/S count a prefix)
    const size_t n = modes.size();
    for (auto* v : { &board.wr, &board.wi, &board.qr, &board.qi, &board.hr, &board.hi, &board.omega, &board.sigma,
                     &board.radLeft, &board.radRight, &board.volume, &board.direct })
        v->assign (n, 0.0f);
    board.phi.assign (n * Z, 0.0f);
    board.phiH.assign (n * Z, 0.0f);
    board.phiSum.assign (n, 0.0f);
    // The solved sideways shapes (the bridge top tilting with the board's slope along the
    // bridge) reach the strings as they are. (They used to be scaled to 0.45, fitted to the
    // recorded tail decays, which was holding a non-passive discretization below its limit.)
    const double sideScale = 1.0;
    board.scratchF.assign (n, 0.0f);
    board.scratchV.assign (n, 0.0f);

    // Sound at the player's seat (or a close mic pair over the strings), worked out with the
    // shapes: both faces of the open board, the short circuit round its edge, the floor and the
    // open lid. One overall scale sets the output units (the same for the frame's modes).
    std::vector<double> rawL (n), rawR (n), rawV (n);
    for (size_t i = 0; i < n; ++i)
    {
        const auto& src = kBoardModes[modes[i].n];
        rawL[i] = settings.closeMics ? src.micLeft : src.radLeft;
        rawR[i] = settings.closeMics ? src.micRight : src.radRight;
        rawV[i] = src.volume;
    }
    // (The seat's octave balance used to be forced to the close pair's below 1 kHz, a fit to the
    // C5 recordings; the seat now hears the solved body as it is.)
    {
        double sum = 0; int count = 0;
        for (size_t i = 0; i < n; ++i)
            if (modes[i].f >= 200.0 && modes[i].f < 800.0) { sum += rawL[i] * rawL[i] + rawR[i] * rawR[i]; count += 2; }
        radScale = count > 0 && sum > 0 ? 0.85 / std::sqrt (sum / count) : 1.0;
        for (size_t i = 0; i < n; ++i) { rawL[i] *= radScale; rawR[i] *= radScale; }
    }
    // The sound leaking out under the board (see buildBody) follows the net volume a mode pushes
    {
        double sum = 0; int count = 0;
        for (size_t i = 0; i < n; ++i)
            if (modes[i].f < 400.0) { sum += rawV[i] * rawV[i]; ++count; }
        const double scale = count > 0 && sum > 0 ? 0.85 / std::sqrt (sum / count) : 0.0;
        for (auto& v : rawV) v *= scale;
    }

    for (size_t i = 0; i < n; ++i)
    {
        const double f = modes[i].f;
        const double w = 2.0 * M_PI * f;
        const double sigma = M_PI * f * modes[i].loss;
        // Unit modal mass: the shapes carry the physical admittance. The mode is discretized
        // with the trapezoidal rule (prewarped so its resonance lands exactly), not by impulse
        // invariance: the trapezoidal admittance is positive real at every frequency, like the
        // physical one, so the strings can be coupled to it at full strength without the
        // coupling feeding energy in. In partial fractions it is a direct term plus a complex
        // one-pole pair: velocity = direct * force + Re[(omega + j sigma) w], w <- q w + h force
        const double w0sq = w * w + sigma * sigma;
        const double K = w / std::tan (0.5 * w * Ts);
        const double a0 = K * K + 2.0 * sigma * K + w0sq, a1 = 2.0 * (w0sq - K * K), a2 = K * K - 2.0 * sigma * K + w0sq;
        const double d = -K / a2, c0 = K - d * a0, c1 = -d * a1;
        const double qr = -a1 / (2.0 * a0), qi = std::sqrt (std::max (0.0, a2 / a0 - qr * qr));
        const double ar = c0 / (2.0 * a0), ai = (-c1 / (2.0 * a0) - ar * qr) / qi;
        const cd pole (qr, qi);
        const cd gain = 2.0 * cd (ar, ai) / cd (w, sigma);
        board.omega[i] = (float) w;
        board.sigma[i] = (float) sigma;
        board.direct[i] = (float) d;
        board.qr[i] = (float) pole.real();
        board.qi[i] = (float) pole.imag();
        board.hr[i] = (float) gain.real();
        board.hi[i] = (float) gain.imag();
        // The bridge's own mass (see buildBody) also shields the strings from board modes
        // above its corner, most in the heavy bass and tenor
        for (int z = 0; z < Z; ++z)
        {
            const double fc = 2000.0 * std::pow (1.22, layoutPosition (z)), ratio = modes[i].f / fc;
            board.phi[(size_t) z * n + i] = (float) (phiScale * zoneCoupling[(size_t) z] * modes[i].phi[(size_t) z] / std::sqrt (1.0 + ratio * ratio));
            board.phiH[(size_t) z * n + i] = (float) (sideScale * phiScale * zoneCoupling[(size_t) z] * kBoardModes[modes[i].n].phiH[z] / std::sqrt (1.0 + ratio * ratio));
            board.phiSum[i] += board.phi[(size_t) z * n + i];
        }
        board.radLeft[i] = (float) rawL[i];
        board.radRight[i] = (float) rawR[i];
        board.volume[i] = (float) rawV[i];
    }

    // Above the solved modes the board is a many-mode system: its modes overlap so much that
    // the bridge sees a smooth, resistive load. Its mean value follows from the modes just
    // below the top (the average real admittance of overlapping modes is pi/2 times their
    // squared shapes per unit of angular frequency), and it carries on above them.
    {
        const double f0 = 0.7 * fTop;
        for (int z = 0; z < Z; ++z)
        {
            double sum = 0;
            for (size_t i = 0; i < n; ++i)
                if (modes[i].f >= f0 && modes[i].f < fTop)
                {
                    const double p = board.phi[(size_t) z * n + i];
                    sum += p * p;
                }
            highAdmittance[(size_t) z] = 0.5 * M_PI * sum / (2.0 * M_PI * (fTop - f0));
            highPassState[(size_t) z] = highPassState2[(size_t) z] = 0;
        }
        highPassCoeff = std::exp (-2.0 * M_PI * fTop / sampleRate);
        // Its sound: overlapping modes add in power, so the sound per unit of bridge velocity
        // is the root of (their radiation squared) over (their shapes squared at the bridge)
        double rad2 = 0, phi2 = 0;
        for (size_t i = 0; i < n; ++i)
            if (modes[i].f >= f0 && modes[i].f < fTop)
            {
                rad2 += 0.5 * (rawL[i] * rawL[i] + rawR[i] * rawR[i]);
                for (int z = 0; z < Z; ++z)
                    phi2 += (double) board.phi[(size_t) z * n + i] * board.phi[(size_t) z * n + i] / Z;
            }
        highRadiation = phi2 > 0 ? std::sqrt (rad2 / phi2) : 0.0;
        highLeftState = highRightState = highLeftState2 = highRightState2 = 0;
    }

    // Overall level: a mezzo-forte middle C peaks around -26 dBFS, fortissimo around -16, so
    // that a full fortissimo chord with the pedal down stays under the plugin's safety clipper
    // (it rounds off peaks above -2 dBFS; a Nocturne render at the old gain pushed 6000
    // samples into it, which was heard as distortion on the loud notes)
    outputGain = 21.5 * 0.74 * 0.5;
    buildBridgeAdmittance();
    buildBody();
}

void ModalPiano::buildBridgeAdmittance()
{
    // The board's motion at each zone in the same sample as a unit force at each zone (vertical
    // zones first, then sideways): the modes' direct feedthrough, summed over the shapes. The
    // many-mode admittance above the solved modes adds its instantaneous part on the vertical
    // diagonal, and the leaning bridge pins (side bearing) fold in on both sides, as the strings
    // see them: the sideways pull pushes up by sideBearing, and the up-and-down motion moves the
    // sideways strings by the same share
    const int n = (int) board.omega.size(), Z = kZones, P = 2 * kZones;
    std::vector<double> y ((size_t) P * P, 0.0);
    std::vector<double> gain ((size_t) n);
    for (int i = 0; i < n; ++i)
        gain[(size_t) i] = (double) board.omega[(size_t) i] * board.hr[(size_t) i] - (double) board.sigma[(size_t) i] * board.hi[(size_t) i] + board.direct[(size_t) i];
    auto row = [&] (int p) { return p < Z ? board.phi.data() + (size_t) p * (size_t) n : board.phiH.data() + (size_t) (p - Z) * (size_t) n; };
    for (int p = 0; p < P; ++p)
        for (int q = p; q < P; ++q)
        {
            const float* a = row (p);
            const float* b = row (q);
            double sum = 0;
            for (int i = 0; i < n; ++i)
                sum += (double) a[i] * b[i] * gain[(size_t) i];
            y[(size_t) p * P + q] = y[(size_t) q * P + p] = sum;
        }
    const double c = highPassCoeff;
    for (int z = 0; z < Z; ++z)
        y[(size_t) z * P + z] += highAdmittance[(size_t) z] * c * c;
    // Y_eff = B^T Y B with B = [[I, sb I], [0, I]] (string forces to board forces)
    bridgeY.assign ((size_t) P * P, 0.0);
    const double sb = sideBearing;
    for (int p = 0; p < P; ++p)
        for (int q = 0; q < P; ++q)
        {
            // (Y B)_pq
            double yb = y[(size_t) p * P + q];
            if (q >= Z)
                yb += sb * y[(size_t) p * P + (q - Z)];
            // (B^T (Y B))_pq: row p of B^T is row p of B transposed
            double v = yb;
            if (p >= Z)
            {
                double yb2 = y[(size_t) (p - Z) * P + q];
                if (q >= Z)
                    yb2 += sb * y[(size_t) (p - Z) * P + (q - Z)];
                v += sb * yb2;
            }
            bridgeY[(size_t) p * P + q] = v;
        }
}

void ModalPiano::buildFrame()
{
    // The cast-iron frame holds the strings' far ends (agraffes, capo bar). Its modes come from
    // the same body model as the board's (FrameModes.inc: the body modes that mostly move the
    // frame, solved with the rim, board and air), with their shapes where each key zone's
    // strings end, their own losses, and the sound they make through the board.
    static_assert (sizeof (kFrameModes[0].phi) / sizeof (float) == kZones, "regenerate FrameModes.inc for this many zones");
        std::vector<int> use;
    for (int i = 0; i < kFrameModeCount; ++i)
        if (kFrameModes[i].f < 0.45 * sampleRate)
            use.push_back (i);
    const size_t n = use.size();
    for (auto* v : { &frame.wr, &frame.wi, &frame.qr, &frame.qi, &frame.hr, &frame.hi, &frame.omega, &frame.sigma, &frame.rad })
        v->assign (n, 0.0f);
    frame.phi.assign (n * kZones, 0.0f);
    for (size_t k = 0; k < n; ++k)
    {
        const auto& src = kFrameModes[use[k]];
        // At least the loss of the frame's bolted seats on the rim (Q of about 100 with its felt pads)
        const double f = src.f, w = 2.0 * M_PI * f, sigma = M_PI * f * std::max (0.01, (double) src.loss);
        const cd pole = std::exp (cd (-sigma, w) * Ts);
        const cd gain = pole * (Ts / w);
        frame.omega[k] = (float) w; frame.sigma[k] = (float) sigma;
        frame.qr[k] = (float) pole.real(); frame.qi[k] = (float) pole.imag();
        frame.hr[k] = (float) gain.real(); frame.hi[k] = (float) gain.imag();
        // Fitted down like the board (see buildSoundboard): at full strength the frame rang as a chord on every key
        for (int z = 0; z < kZones; ++z)
            frame.phi[k * kZones + (size_t) z] = (float) (0.2 * src.phi[z]);
        // The bars are narrow, so they push little air themselves: most of their sound goes out through the board
        frame.rad[k] = (float) (0.3 * radScale * src.rad);
    }
    frame.accel.fill (0); frame.lastVel.fill (0); frame.force.fill (0);
    frame.knock = 0;
    frame.knockDecay = std::exp (-1.0 / (0.0004 * sampleRate));
}

void ModalPiano::runFrame (double& left, double& right)
{
    // Far-end force of every sounding string: T dy/dx at x = 0 is T n pi / L per unit of mode
    // n, the bridge-end coupling without its alternating sign. The frame's motion where the
    // string ends moves that end back, which the modes feel as beta (-1)^(n+1) times its
    // acceleration.
    frame.force.fill (0);
    for (int index : frameList)
    {
        const auto& st = strings[(size_t) index];
        const float a = (float) (frame.accel[(size_t) st.zone] * settings.farEnd);
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
        frame.force[(size_t) st.zone] += f * settings.farEnd;
    }
    // The hammer's blow also travels through the strings' far ends into the iron plate: a
    // short, sharp push that sets the frame ringing for a moment (the metallic shock of a
    // hard treble note), spread along the front bar
    for (auto& fz : frame.force)
        fz += frame.knock / kZones;
    frame.knock *= frame.knockDecay;

    std::array<double, kZones> v {};
    double out = 0;
    const size_t n = frame.omega.size();
    for (size_t i = 0; i < n; ++i)
    {
        const float* ph = frame.phi.data() + i * kZones;
        double F = 0;
        for (int z = 0; z < kZones; ++z)
            F += ph[z] * frame.force[(size_t) z];
        const float x = frame.wr[i], y = frame.wi[i];
        const float nr = frame.qr[i] * x - frame.qi[i] * y + frame.hr[i] * (float) F;
        const float ni = frame.qr[i] * y + frame.qi[i] * x + frame.hi[i] * (float) F;
        frame.wr[i] = nr; frame.wi[i] = ni;
        const float vel = frame.omega[i] * nr - frame.sigma[i] * ni;
        for (int z = 0; z < kZones; ++z)
            v[(size_t) z] += ph[z] * vel;
        out += frame.rad[i] * vel;
    }
    for (int z = 0; z < kZones; ++z)
    {
        // The far end's own give: the agraffes in the bass hold the string firmly, the capo
        // bar in the treble and the bearing on the pinblock felt let it move and lose more
        const double farR = 3e-6 * (0.7 + 0.6 * z / (kZones - 1.0));
        const double vz = v[(size_t) z] + farR * frame.force[(size_t) z];
        frame.accel[(size_t) z] = (vz - frame.lastVel[(size_t) z]) * sampleRate;
        frame.lastVel[(size_t) z] = vz;
    }
    left += out;
    right += 0.9 * out;
}

void ModalPiano::buildBody()
{
    // Air under the board: its resonances, its push back on the board and the sound leaking
    // out of the open sides and off the floor are all in the body solve (BoardModes.inc), so
    // there is no separate leak path here any more: it counted the rocking modes twice.
    airModes.clear();
    auto addAir = [&] (double f, double q, double g)
    {
        const double w = 2.0 * M_PI * f / sampleRate;
        const double alpha = std::sin (w) / (2.0 * q);
        const double a0 = 1.0 + alpha;
        airModes.push_back ({ (float) (alpha / a0), (float) (-2.0 * std::cos (w) / a0), (float) ((1.0 - alpha) / a0),
                              0, 0, 0, 0, (float) g });
    };
    (void) addAir;
    airPressure = 0;

    // The lid's reflection is in the body tables: each mode's sound at the seat and the mics
    // includes the mirror image off the finite lid (build_board.py), so nothing is added here.

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
        // Many overlapping modes: a dense, random ripple of a few dB, different for each ear
        for (double f = 5200.0; f < 0.42 * sampleRate; f *= 1.12)
        {
            highBoardLeft.push_back (peaking (f * (1.0 + 0.04 * rng.symmetric()), 3.0 * rng.symmetric(), 4.0));
            highBoardRight.push_back (peaking (f * (1.0 + 0.04 * rng.symmetric()), 3.0 * rng.symmetric(), 4.0));
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

    // The lid is in the body tables: build_board.py solves it as a hinged, propped plate over
    // the case air, coupled to the board and radiating to the listener with the other modes.

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
    zoneSideForce.fill (0); zoneSideAccel.fill (0); lastZoneSideVel.fill (0); zoneSideBoardVel.fill (0);
    lastHorizontalVelocity = horizontalAcceleration = 0;
    resistiveLpLeft = resistiveLpRight = resistiveLp2Left = resistiveLp2Right = resistiveLp3Left = resistiveLp3Right = 0;
    buildBody();
    buildFrame();
    dcLeft = dcRight = dcInLeft = dcInRight = 0;
    sustain = false;
    pedalLevel = 0;
    pedalGrip = 1.0;
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

std::complex<double> ModalPiano::restingPole (const String& s, const KeyState& ks, size_t j) const
{
    if (! s.damped)
        return cd (m.prFree[j], m.piFree[j]);
    // The felt may be settling after landing, or held partly clear by a half-pressed pedal:
    // interpolate the poles in the log domain, so decay and frequency move together from the
    // free values to the damped ones
    const double settle = ks.damperRamp > 0 ? 1.0 - ks.damperRamp / (kDamperRampMs * 0.001 * sampleRate) : 1.0;
    const double g = settle * pedalGrip;
    if (g >= 1.0)
        return cd (m.prDamp[j], m.piDamp[j]);
    if (g <= 0.0)
        return cd (m.prFree[j], m.piFree[j]);
    const cd lf = std::log (cd (m.prFree[j], m.piFree[j]));
    const cd ld = std::log (cd (m.prDamp[j], m.piDamp[j]));
    return std::exp (lf + g * (ld - lf));
}

void ModalPiano::updateDamping (int keyIndex)
{
    const auto& key = keys[(size_t) keyIndex];
    auto& ks = keyState[(size_t) keyIndex];
    const bool damped = key.hasDamper && ! ks.down && ks.damperPending <= 0 && pedalGrip > 0.0 && ! ks.sostenuto;

    forEachString (ks, [&] (String& s)
    {
        if (s.damped != damped || s.activity == Activity::Off)
        {
            s.damped = damped;
            forModes (s, s.count, [&] (size_t j)
            {
                const cd pole = restingPole (s, ks, j);
                m.pr[j] = (float) pole.real();
                m.pi[j] = (float) pole.imag();
            });
        }
        // Undamped strings listen to the bridge
        if (! damped && s.activity == Activity::Off)
            setActivity (s, Activity::Sympathetic);
    });
}

void ModalPiano::setKeyDown (int midiNote, bool down, double releaseSpeed)
{
    const int k = midiNote - kLowestKey;
    if (k < 0 || k >= (int) keys.size())
        return;
    auto& ks = keyState[(size_t) k];
    if (! down && ks.down)
    {
        // The key takes a moment to rise far enough for the damper to reach the strings:
        // 60 ms when it is eased up, less the faster the finger leaves. But the damper only
        // falls under its own weight (the key's back end drops away from it, it is not pushed),
        // so from its 5 mm lift it cannot land in less than sqrt(2 h / g) = 32 ms however the
        // key is let go. Re-press it before then (a fast repeat or a half-stroke) and the
        // damper never lands; a pedal that comes down before then catches nothing.
        ks.releaseSpeed = (float) std::clamp (releaseSpeed, 0.0, 1.0);
        constexpr double kDamperDrop = 0.032;
        const double fall = std::max (kDamperDrop, 0.06 - 0.11 * ks.releaseSpeed + 0.06 * ks.releaseSpeed * ks.releaseSpeed);
        ks.releasedAt = sampleClock;
        ks.damperPending = std::max (1, (int) (fall * sampleRate));
        ks.down = false;
        return;
    }
    ks.down = down;
    ks.damperPending = 0;
    updateDamping (k);
    checkActivity();
}

void ModalPiano::setSustainLevel (double level)
{
    level = std::clamp (level, 0.0, 1.0);
    // The pedal's first quarter of travel takes up the lost motion; from there the damper
    // tray lifts the felts off gradually, and they are clear of the strings by about 0.65.
    // The felt's grip falls away faster than its pressure (squared), so the half-pedal
    // range (CC64 50-80) runs from a quick damping to a barely slowed ring.
    const double lift = std::clamp ((0.65 - level) / 0.4, 0.0, 1.0);
    const double grip = lift * lift;
    const bool edge = (grip <= 0.0) != (pedalGrip <= 0.0) || (grip >= 1.0) != (pedalGrip >= 1.0);
    if (! edge && std::abs (grip - pedalGrip) < 0.02)
        return;
    pedalLevel = level;
    pedalGrip = grip;
    sustain = grip <= 0.0;
    for (size_t k = 0; k < keys.size(); ++k)
        updateDamping ((int) k);
    // Half pedalling moves the grip of every damper that is already down
    if (grip > 0.0)
        for (size_t k = 0; k < keys.size(); ++k)
            forEachString (keyState[k], [&] (String& s)
            {
                if (! s.damped || s.activity == Activity::Off)
                    return;
                forModes (s, s.count, [&] (size_t j)
                {
                    const cd pole = restingPole (s, keyState[k], j);
                    m.pr[j] = (float) pole.real();
                    m.pi[j] = (float) pole.imag();
                });
            });
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

    ks.hammer.softFelt = softPedal;
    // Felt memory between blows: felt squeezed by a blow takes a fraction of a second to
    // spring back, so a key repeated quickly meets the strings with thinner, harder felt.
    // The harder the earlier blows, the more is left
    {
        const double since = ks.struckAt < 0 ? 1e9 : (sampleClock - ks.struckAt) / sampleRate;
        ks.feltSet = ks.feltSet * std::exp (-since / 0.2);
        ks.hammer.feltSet = ks.feltSet;
        ks.feltSet = std::min (1.0, ks.feltSet + 0.5 * std::clamp (hammerVelocity / 5.0, 0.0, 1.0));
        ks.struckAt = sampleClock;
    }

    ks.hammer.active = true;
    ks.hammer.position = nearest - letOff; // just let off, flying free
    ks.hammer.velocity = std::max (0.05, hammerVelocity);
    ks.hammer.shankPos = nearest - letOff;
    ks.hammer.shankVel = ks.hammer.velocity;
    ks.hammer.samplesSinceStrike = 0;
    ks.hammer.samplesApart = 0;
    for (int i = 0; i < ks.numStrings; ++i)
    {
        const auto& s = strings[(size_t) (ks.firstString + i)];
        double u = 0;
        forModes (s, s.count, [&] (size_t j) { u += m.phi[j] * m.zi[j]; });
        ks.hammer.lastCompression[(size_t) i] = nearest - letOff - u - ks.hammer.gap[(size_t) i];
        ks.hammer.forceSum[(size_t) i] = 0;
        ks.hammer.lastFelt[(size_t) i] = 0;
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
        const bool full = s.activity == Activity::Full;
        segments.push_back ({ s.first, s.sympatheticCount, s.panLeft, s.panRight, side, s.zone, full ? 1e9f : s.sympTopF });
        if (full && s.count > s.sympatheticCount)
            segments.push_back ({ s.highFirst, s.count - s.sympatheticCount, s.panLeft, s.panRight, side, s.zone, 1e9f });
    }
    std::sort (segments.begin(), segments.end(), [] (const Segment& a, const Segment& b) { return a.first < b.first; });

    zoneLive.fill (0);
    zoneSideLive.fill (0);
    std::array<float, kZones> topV {}, topS {};
    for (const auto& g : segments)
    {
        zoneLive[(size_t) g.zone] = 1;
        if (g.side)
            zoneSideLive[(size_t) g.zone] = 1;
        auto& top = g.side ? topS[(size_t) g.zone] : topV[(size_t) g.zone];
        top = std::max (top, g.topF);
    }
    // A string under the hammer is stepped by the hammer code, not from a segment: its zone
    // still needs the whole board
    for (int k : hammerList)
    {
        zoneLive[(size_t) keyZone[(size_t) k]] = 1;
        topV[(size_t) keyZone[(size_t) k]] = 1e9f;
    }
    for (int z = 0; z < kZones; ++z)
    {
        // modes up to a third above the highest resonating partial (the board modes just above
        // it still shape the bridge there); all of them under a struck string
        auto count = [&] (float topF) {
            const float wTop = (float) (2.0 * M_PI) * topF * 1.33f;
            return (int) (std::upper_bound (board.omega.begin(), board.omega.end(), wTop) - board.omega.begin());
        };
        zoneModesV[(size_t) z] = zoneLive[(size_t) z] ? count (topV[(size_t) z]) : 0;
        zoneModesS[(size_t) z] = zoneSideLive[(size_t) z] ? count (topS[(size_t) z]) : 0;
    }
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
                      * (1.0 - 0.35 * h.softFelt) * (1.0 + 0.4 * h.feltSet);
    const double p = key.feltExponent;   // the unpacked felt beside the grooves is softer (q0 above), not a different law
    const double coreDepth = key.outerFeltDepth * 1000.0 * settings.outerFeltScale * (1.0 - 0.6 * ageFactor) * (1.0 - 0.5 * groove) * (1.0 - 0.3 * h.feltSet); // [mm]
    // Treble hammers are smaller and their outer felt firmer, or soft notes there would vanish
    const double outer = settings.outerFelt * std::pow (10.0, std::max (0.0, (keyIndex - 39) / 48.0) * trebleFirmness);
    // The felt's loss (see ScaleDesign): the viscous term is the elastic one's rate of change
    // times the key's felt time. It stiffens the felt as it is squeezed and softens it as it
    // lets go, so the hammer gets back less than it gave; the felt never pulls on the string.
    const double tau = key.feltDamping;
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
            const double viscous = tau * (elastic - h.lastFelt[(size_t) i]) / dt;
            h.lastFelt[(size_t) i] = elastic;
            const double f = compression > 0 ? q0 * h.feltScale[(size_t) i] * std::max (0.0, elastic + viscous) : 0.0;

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

        // Hammer flight (symplectic Euler). The head (Chabassier's measured head masses, 11 g
        // at A0 to 5.8 g at C8) rides on the end of the shank, a hornbeam cantilever (E 13 GPa,
        // 6.5 mm thick, 125 mm from the flange, about 3.3 g, of which a third counts at its
        // end as it swings about the flange). Its bending stiffness 3EI/L^3 is about 1.75 kN/m,
        // so with the head on its end it bends at 60-100 Hz: during the 1-2 ms blow the hammer
        // is rigid for the strings (contact durations as published: about 2 ms at C4, 1.2 ms
        // at C6), and the bending only shapes how the hammer leaves (Askenfelt & Jansson 1991).
        if (settings.shankFlex > 0)
        {
            const double mh = key.hammerMass, ms = 1.1e-3;
            const double EI = 13e9 * M_PI * std::pow (6.5e-3, 4) / 64.0;
            const double kb = settings.shankFlex * 3.0 * EI / std::pow (0.125, 3);
            const double cb = 2.0 * 0.1 * std::sqrt (kb * mh * ms / (mh + ms));
            const double fBend = kb * (h.shankPos - h.position) + cb * (h.shankVel - h.velocity);
            h.velocity += (fBend - hammerForce) / mh * dt;
            h.shankVel += (-fBend) / ms * dt;
            h.position += h.velocity * dt;
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
            forModes (s, s.count, [&] (size_t j)
            {
                const cd base = restingPole (s, ks, j);
                const float pr0 = (float) base.real(), pi0 = (float) base.imag();
                const float w = m.omega[j];
                const float y = w * shift;
                const float a = 1.0f - w * decay - 0.5f * y * y;
                m.pr[j] = a * pr0 - y * pi0;
                m.pi[j] = a * pi0 + y * pr0;
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
            if (keys[k].hasDamper && pedalGrip > 0.0 && ! ks.sostenuto)
            {
                // The felt lands harder from a snapped-off key than from one eased up
                queueAction (2, 0, 0.3 * (0.4 + 2.4 * ks.releaseSpeed * ks.releaseSpeed) * pedalGrip, (int) k);
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
                ks.damperRamp = (int) (kDamperRampMs * 0.001 * sampleRate);
            }
            updateDamping ((int) k);
            checkActivity();
        }
    }
    // A landed damper's felt settles onto the strings over a few tens of milliseconds, so a
    // quick re-pedal after a release (a 30 ms gap is common) only grazes the note instead of
    // cutting it by 10 dB. Its grip ramps up from nothing to the full damped poles.
    for (size_t k = 0; k < keyState.size(); ++k)
    {
        auto& ks = keyState[k];
        if (ks.damperRamp <= 0)
            continue;
        ks.damperRamp = std::max (0, ks.damperRamp - numSamples);
        forEachString (ks, [&] (String& s)
        {
            if (! s.damped) { ks.damperRamp = 0; return; }
            forModes (s, s.count, [&] (size_t j)
            {
                const cd pole = restingPole (s, ks, j);
                m.pr[j] = (float) pole.real();
                m.pi[j] = (float) pole.imag();
            });
        });
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
        zoneForce.fill (0);
        zoneSideForce.fill (0);
        double forceSum = 0, forceLeft = 0, forceRight = 0;  // vertical
        double sideSum = 0, sideLeft = 0, sideRight = 0;     // horizontal

        // Strings without a hammer on them: one step per sample. The bridge acceleration they
        // feel this sample is not known yet (it depends on their own force on the bridge), so each
        // run folds last sample's acceleration into its state first and reports its force without
        // this sample's input (free) and its gain on that input; the bridge solve below closes
        // the loop in the same sample
        // (the runs can be laid out afresh later in this sample, by a hammer leaving or a string
        // falling quiet, so what each run needs for the solve is kept here, not in the run)
        const size_t stepped = runs.size();
        if (runGain.size() < stepped) { runGain.resize (stepped); runPanLeft.resize (stepped); runPanRight.resize (stepped); runDim.resize (stepped); }
        std::array<double, 2 * kZones> implicitGain {};
        for (size_t ri = 0; ri < stepped; ++ri)
        {
            const auto& r = runs[ri];
            const int j = r.first;
            float g = 0;
            const float f = stepModesDeferred (zr + j, zi + j, pr + j, pi + j, bgr + j, bgi + j, fk + j, r.count,
                                               (float) (r.side ? zoneSideAccel[(size_t) r.zone] : zoneAcceleration[(size_t) r.zone]), g);
            runGain[ri] = g;
            runPanLeft[ri] = r.panLeft;
            runPanRight[ri] = r.panRight;
            runDim[ri] = r.side ? kZones + r.zone : r.zone;
            implicitGain[(size_t) runDim[ri]] += g;
            if (r.side)
            {
                sideSum += f;
                zoneSideForce[(size_t) r.zone] += f;
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
            actionSpread = (forceSum - before) / kZones;
        }

        // Soundboard and frame modes. Each mode is pushed by every zone's force times its
        // shape there, and moves every zone by the same shape. The free strings' forces depend on
        // the bridge's acceleration this sample, which depends on those forces: the board first
        // moves on its own (pass 1), the live zones' small linear system is solved for the
        // accelerations (the board's same-sample give per unit force, bridgeY, against the
        // strings' gain on their input), and the forces are then applied (pass 2). No sample of
        // delay is left in the loop, so the full sideways lever of a tall bridge stays stable.
        double radiatedLeft = 0, radiatedRight = 0, volumeVelocity = 0;
        zoneBoardVelocity.fill (0);
        {
            // Laid out zone by zone (phi[z * modes + i]) so each loop runs over the modes
            const int nb = (int) boardModes;
            float* __restrict F = board.scratchF.data();
            float* __restrict V = board.scratchV.data();
            const float* __restrict phi = board.phi.data();
            const float* __restrict phH = board.phiH.data();
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
            const float* __restrict dr = board.direct.data();
            // Pass 1: the modes' free motion this sample
            float l = 0, r = 0, vol = 0;
            for (int i = 0; i < nb; ++i)
            {
                const float x = wr[i], y = wi[i];
                const float nr = qr[i] * x - qi[i] * y;
                const float ni = qr[i] * y + qi[i] * x;
                wr[i] = nr;
                wi[i] = ni;
                const float vel = om[i] * nr - sg[i] * ni;
                V[i] = vel;
                l += rl[i] * vel;
                r += rrt[i] * vel;
                vol += bv[i] * vel;
            }
            // The free bridge motion at each zone that has strings being run (and sideways,
            // where sideways strings ring): four shapes per pass over the modes' velocities
            float dotOut[2 * kZones];
            const float* rows[2 * kZones];
            int lens[2 * kZones];
            int nd = 0;
            for (int z = 0; z < kZones; ++z)
            {
                if (zoneLive[(size_t) z]) { rows[nd] = phi + (size_t) z * (size_t) nb; lens[nd++] = zoneModesV[(size_t) z]; }
                // (wherever a sideways string runs, not only where one already pushes: a
                // resting sideways string must hear the board to resonate in sympathy)
                if (zoneSideLive[(size_t) z]) { rows[nd] = phH + (size_t) z * (size_t) nb; lens[nd++] = zoneModesS[(size_t) z]; }
            }
            int d0 = 0;
            for (; d0 + 4 <= nd; d0 += 4)
            {
                const float* __restrict p0 = rows[d0]; const float* __restrict p1 = rows[d0 + 1];
                const float* __restrict p2 = rows[d0 + 2]; const float* __restrict p3 = rows[d0 + 3];
                float c0 = 0, c1 = 0, c2 = 0, c3 = 0;
                const int common = std::min (std::min (lens[d0], lens[d0 + 1]), std::min (lens[d0 + 2], lens[d0 + 3]));
                for (int i = 0; i < common; ++i)
                {
                    const float v = V[i];
                    c0 += p0[i] * v; c1 += p1[i] * v; c2 += p2[i] * v; c3 += p3[i] * v;
                }
                float* cs[4] = { &c0, &c1, &c2, &c3 };
                for (int q = 0; q < 4; ++q)
                {
                    const float* __restrict p = rows[d0 + q];
                    float c = 0;
                    for (int i = common; i < lens[d0 + q]; ++i)
                        c += p[i] * V[i];
                    *cs[q] += c;
                }
                dotOut[d0] = c0; dotOut[d0 + 1] = c1; dotOut[d0 + 2] = c2; dotOut[d0 + 3] = c3;
            }
            for (; d0 < nd; ++d0)
            {
                const float* __restrict p0 = rows[d0];
                float c0 = 0;
                for (int i = 0; i < lens[d0]; ++i)
                    c0 += p0[i] * V[i];
                dotOut[d0] = c0;
            }

            // The bridge this sample, zone by zone: p runs over the live vertical zones then the
            // live sideways ones. The strings give f_p = K_p + G_p a_p (K: every force already
            // known, from hammered strings, the lengthwise modes and the free strings' own
            // motion; G: the free strings' gain on this sample's acceleration); the board answers
            // v_p = vFree_p + sum_q Y_pq f_q, with the many-mode admittance above the solved
            // modes and the leaning pins (side bearing) folded into Y and vFree; and
            // a_p = (v_p - last v_p) * rate. Solved for a.
            int dims[2 * kZones];
            int np = 0;
            double vFree[2 * kZones];
            double highFree[kZones];
            {
                const double c = highPassCoeff;
                int di = 0;
                for (int z = 0; z < kZones; ++z)
                {
                    // The smooth admittance's two-pole high-pass, split into what its memory gives
                    // on its own and (in bridgeY) what this sample's force adds at once
                    highFree[(size_t) z] = highAdmittance[(size_t) z] * (-c * c * highPassState[(size_t) z] - c * highPassState2[(size_t) z]);
                    // (the dot products come zone by zone, up-and-down then sideways, as the rows went in)
                    if (zoneLive[(size_t) z])
                    {
                        vFree[(size_t) z] = dotOut[di++] + highFree[(size_t) z];
                        dims[np++] = z;
                    }
                    if (zoneSideLive[(size_t) z])
                        vFree[(size_t) (kZones + z)] = dotOut[di++] + sideBearing * vFree[(size_t) z];
                }
                for (int z = 0; z < kZones; ++z)
                    if (zoneSideLive[(size_t) z])
                        dims[np++] = kZones + z;
            }
            double accel[2 * kZones] = {};
            {
                // (I - rate Y G) a = rate (vFree + Y K - vLast)
                double A[4 * kZones * kZones], rhs[2 * kZones];
                const double rate = sampleRate;
                for (int u = 0; u < np; ++u)
                {
                    const int p = dims[u];
                    double yk = 0;
                    for (int w = 0; w < np; ++w)
                    {
                        const int q = dims[w];
                        const double y = bridgeY[(size_t) p * 2 * kZones + (size_t) q];
                        yk += y * (q < kZones ? zoneForce[(size_t) q] : zoneSideForce[(size_t) (q - kZones)]);
                        A[u * np + w] = (u == w ? 1.0 : 0.0) - rate * y * implicitGain[(size_t) q];
                    }
                    const double last = p < kZones ? lastZoneVelocity[(size_t) p] : lastZoneSideVel[(size_t) (p - kZones)];
                    rhs[u] = rate * (vFree[(size_t) p] + yk - last);
                }
                solveDense (A, rhs, np);
                for (int u = 0; u < np; ++u)
                    accel[(size_t) dims[u]] = rhs[u];
            }
            // The free strings' forces, complete
            for (size_t ri = 0; ri < stepped; ++ri)
            {
                const int p = runDim[ri];
                const double df = runGain[ri] * accel[(size_t) p];
                if (df == 0.0)
                    continue;
                if (p >= kZones)
                {
                    sideSum += df;
                    zoneSideForce[(size_t) (p - kZones)] += df;
                    sideLeft += df * runPanLeft[ri];
                    sideRight += df * runPanRight[ri];
                }
                else
                {
                    forceSum += df;
                    forceLeft += df * runPanLeft[ri];
                    forceRight += df * runPanRight[ri];
                    zoneForce[(size_t) p] += df;
                }
            }
            // The bridge's motion, as the solve saw it, fed to every string next sample
            for (int u = 0; u < np; ++u)
            {
                const int p = dims[u];
                double v = vFree[(size_t) p];
                for (int w = 0; w < np; ++w)
                {
                    const int q = dims[w];
                    v += bridgeY[(size_t) p * 2 * kZones + (size_t) q] * (q < kZones ? zoneForce[(size_t) q] : zoneSideForce[(size_t) (q - kZones)]);
                }
                if (p < kZones)
                {
                    zoneBoardVelocity[(size_t) p] = v - highFree[(size_t) p];
                    zoneAcceleration[(size_t) p] = accel[(size_t) p];
                    lastZoneVelocity[(size_t) p] = v;
                }
                else
                {
                    zoneSideBoardVel[(size_t) (p - kZones)] = v;
                    zoneSideAccel[(size_t) (p - kZones)] = accel[(size_t) p];
                    lastZoneSideVel[(size_t) (p - kZones)] = v;
                }
            }
            // Side bearing: the strings' sideways pull on the leaning bridge pins has an upward part
            if (sideBearing > 0)
                for (int z = 0; z < kZones; ++z)
                    if (zoneSideForce[(size_t) z] != 0.0)
                        zoneForce[(size_t) z] += sideBearing * zoneSideForce[(size_t) z];
            // Zones without running strings: the many-mode admittance still answers the forces
            // there, and nothing else moves them
            for (int z = 0; z < kZones; ++z)
            {
                const double fz = zoneForce[(size_t) z] + actionSpread;
                auto& hp = highPassState[(size_t) z];
                hp = fz + highPassCoeff * (hp - fz);
                auto& hp2 = highPassState2[(size_t) z];
                const double h1 = fz - hp;
                hp2 = h1 + highPassCoeff * (hp2 - h1);
                if (zoneLive[(size_t) z])
                {
                    if (! zoneSideLive[(size_t) z])
                    {
                        // Sideways: only the up-and-down motion that the leaning pins turn sideways
                        const double vs = sideBearing * lastZoneVelocity[(size_t) z];
                        zoneSideAccel[(size_t) z] = (vs - lastZoneSideVel[(size_t) z]) * sampleRate;
                        lastZoneSideVel[(size_t) z] = vs;
                    }
                    continue;
                }
                const double v = highAdmittance[(size_t) z] * (h1 - hp2);
                zoneAcceleration[(size_t) z] = (v - lastZoneVelocity[(size_t) z]) * sampleRate;
                lastZoneVelocity[(size_t) z] = v;
                zoneSideAccel[(size_t) z] = -lastZoneSideVel[(size_t) z] * sampleRate;
                lastZoneSideVel[(size_t) z] = 0;
            }

            // Pass 2: every zone's push (and sideways pull) times its shape, four rows per pass
            // so F is read and written a quarter as often
            {
                // The action's knocks reach the board through the frame, evenly along the bridge
                const float* __restrict ps = board.phiSum.data();
                const float sp = (float) actionSpread;
                for (int i = 0; i < nb; ++i)
                    F[i] = ps[i] * sp;
            }
            float amounts[2 * kZones];
            int nr = 0;
            for (int z = 0; z < kZones; ++z)
            {
                const float fz = (float) zoneForce[(size_t) z];
                if (fz != 0.0f) { rows[nr] = phi + (size_t) z * (size_t) nb; amounts[nr] = fz; lens[nr++] = zoneLive[(size_t) z] ? zoneModesV[(size_t) z] : nb; }
                // The strings' sideways pull tilts the bridge (its shapes along the bridge)
                const float sz = (float) zoneSideForce[(size_t) z];
                if (sz != 0.0f) { rows[nr] = phH + (size_t) z * (size_t) nb; amounts[nr] = sz; lens[nr++] = zoneModesS[(size_t) z]; }
            }
            int r0 = 0;
            for (; r0 + 4 <= nr; r0 += 4)
            {
                const float* __restrict p0 = rows[r0]; const float* __restrict p1 = rows[r0 + 1];
                const float* __restrict p2 = rows[r0 + 2]; const float* __restrict p3 = rows[r0 + 3];
                const float a0 = amounts[r0], a1 = amounts[r0 + 1], a2 = amounts[r0 + 2], a3 = amounts[r0 + 3];
                const int common = std::min (std::min (lens[r0], lens[r0 + 1]), std::min (lens[r0 + 2], lens[r0 + 3]));
                for (int i = 0; i < common; ++i)
                    F[i] += p0[i] * a0 + p1[i] * a1 + p2[i] * a2 + p3[i] * a3;
                for (int q = 0; q < 4; ++q)
                {
                    const float* __restrict p = rows[r0 + q];
                    const float a = amounts[r0 + q];
                    for (int i = common; i < lens[r0 + q]; ++i)
                        F[i] += p[i] * a;
                }
            }
            for (; r0 < nr; ++r0)
            {
                const float* __restrict p0 = rows[r0];
                const float a0 = amounts[r0];
                for (int i = 0; i < lens[r0]; ++i)
                    F[i] += p0[i] * a0;
            }
            for (int i = 0; i < nb; ++i)
            {
                const float f = F[i];
                wr[i] += hr[i] * f;
                wi[i] += hi[i] * f;
                const float dv = (om[i] * hr[i] - sg[i] * hi[i] + dr[i]) * f;
                l += rl[i] * dv;
                r += rrt[i] * dv;
                vol += bv[i] * dv;
            }
            radiatedLeft = l;
            radiatedRight = r;
            volumeVelocity = vol;
        }

        // Sound: soundboard modes plus the direct (dense high-frequency) bridge path
        // (the many-mode board above the solved modes, then the sideways motion)
        double meanHigh = 0;
        for (int z = 0; z < kZones; ++z) meanHigh += highAdmittance[(size_t) z];
        meanHigh /= kZones;
        highLeftState = forceLeft + highPassCoeff * (highLeftState - forceLeft);
        highRightState = forceRight + highPassCoeff * (highRightState - forceRight);
        const double hl = forceLeft - highLeftState, hr = forceRight - highRightState;
        highLeftState2 = hl + highPassCoeff * (highLeftState2 - hl);
        highRightState2 = hr + highPassCoeff * (highRightState2 - hr);
        const double directLeft = highRadiation * meanHigh * (hl - highLeftState2);
        const double directRight = highRadiation * meanHigh * (hr - highRightState2);
        float dl = (float) directLeft, dr = (float) directRight;
        for (auto& d : diffLeft) dl = d.process (dl);
        for (auto& d : diffRight) dr = d.process (dr);
        for (auto& b : highBoardLeft) dl = b.process (dl);
        for (auto& b : highBoardRight) dr = b.process (dr);
        resistiveLpLeft = dl + lpCoeff * (resistiveLpLeft - dl);
        resistiveLpRight = dr + lpCoeff * (resistiveLpRight - dr);
        // Two more poles: the many-mode board's sound falls off steeply above ~8 kHz (the
        // recordings have ~20 dB less there than one pole leaves, heard as a scuff on loud treble notes)
        resistiveLp2Left = resistiveLpLeft + lp2Coeff * (resistiveLp2Left - resistiveLpLeft);
        resistiveLp2Right = resistiveLpRight + lp2Coeff * (resistiveLp2Right - resistiveLpRight);
        resistiveLp3Left = resistiveLp2Left + lp2Coeff * (resistiveLp3Left - resistiveLp2Left);
        resistiveLp3Right = resistiveLp2Right + lp2Coeff * (resistiveLp3Right - resistiveLp2Right);
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

        double outL = radiatedLeft + resistiveLp3Left + air + frameLeft;
        double outR = radiatedRight + resistiveLp3Right + 0.9 * air + frameRight;

        // Directivity: the high part beams out past the open side
        beamLpLeft = outL + beamCoeff * (beamLpLeft - outL);
        beamLpRight = outR + beamCoeff * (beamLpRight - outR);
        outL -= 0.2 * (outL - beamLpLeft);
        outR += 0.25 * (outR - beamLpRight);

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
