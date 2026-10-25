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

    // Deterministic pseudo-random numbers in [0, 1)
    struct Rng
    {
        uint32_t s;
        explicit Rng (uint32_t seed) : s (seed * 747796405u + 2891336453u) {}
        double next()
        {
            s ^= s << 13; s ^= s >> 17; s ^= s << 5;
            return (s & 0xffffff) / double (0x1000000);
        }
    };

    // Internal (string material and air) losses, as a decay rate in 1/s for a partial
    // at frequency f, on a string whose fundamental is f1. Lower strings ring longer, and
    // high partials die faster. Bridge losses come on top, from the coupling itself.
    double internalDecay (double f, double f1, bool wound)
    {
        const double octavesAboveA0 = std::log2 (std::max (f1, 27.5) / 27.5);
        const double b1 = 0.15 + 0.35 * octavesAboveA0 / 7.25;
        const double b2 = wound ? 2.5e-7 : 1.2e-7;
        return b1 + b2 * f * f;
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

    // Extra decay with the felt damper resting on the string
    double damperDecay (double f, double strength)
    {
        return strength * 10.0 * (1.0 + f / 2000.0);
    }
}

void ModalPiano::Modes::resize (size_t n)
{
    for (auto* v : { &zr, &zi, &pr, &pi, &prFree, &piFree, &prDamp, &piDamp, &bgr, &bgi, &fk, &beta, &phi,
                     &prOs, &piOs, &gOsr, &gOsi, &sigmaInternal, &omega })
        v->assign (n, 0.0f);
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

    // Output filters: gentle top-end roll-off on the direct bridge path, and a DC blocker
    lpCoeff = std::exp (-2.0 * M_PI * 7000.0 / sampleRate);
    dcCoeff = std::exp (-2.0 * M_PI * 20.0 / sampleRate);

    fullList.reserve (strings.size());
    hammerList.reserve (keys.size());
    reset();
}

void ModalPiano::buildModes()
{
    strings.clear();
    keyState.assign (keys.size(), {});

    const double fMax = std::min (0.45 * sampleRate, 16000.0);
    const int os = std::max (1, settings.oversampling);

    // First pass: count modes
    struct Plan { int key, s, count; };
    std::vector<Plan> plans;
    size_t total = 0;
    for (size_t k = 0; k < keys.size(); ++k)
    {
        keyState[k].firstString = (int) plans.size();
        keyState[k].numStrings = keys[k].numStrings;
        for (int s = 0; s < keys[k].numStrings; ++s)
        {
            const auto& str = keys[k].strings[(size_t) s];
            const double fFlex = str.fundamental / std::sqrt (1.0 + str.inharmonicity);
            const int kMax = str.wound ? 110 : 70;
            int count = 0;
            for (int n = 1; n <= kMax; ++n)
            {
                if (n * fFlex * std::sqrt (1.0 + str.inharmonicity * n * n) > fMax)
                    break;
                count = n;
            }
            count = std::max (count, 1);
            plans.push_back ({ (int) k, s, count });
            total += (size_t) count;
        }
    }

    m.resize (total);
    strings.resize (plans.size());

    // Second pass: fill in each mode
    int offset = 0;
    for (size_t i = 0; i < plans.size(); ++i)
    {
        const auto& plan = plans[i];
        const auto& key = keys[(size_t) plan.key];
        const auto& str = key.strings[(size_t) plan.s];
        auto& s = strings[i];

        s.key = plan.key;
        s.first = offset;
        s.count = plan.count;
        s.sympatheticCount = std::min (plan.count, key.hasDamper ? 12 : 8);
        s.mass = (float) (str.linearDensity * str.length / 2.0);

        // Pan with the keyboard, as heard from the player's seat
        const double position = (key.midiNote - kLowestKey) / double (kHighestKey - kLowestKey);
        const double angle = M_PI / 4 + (position - 0.5) * 0.9;
        s.panLeft = (float) std::cos (angle);
        s.panRight = (float) std::sin (angle);

        const double fFlex = str.fundamental / std::sqrt (1.0 + str.inharmonicity);
        const double L = str.length;
        const double strike = key.strikeRatio;
        for (int n = 1; n <= plan.count; ++n)
        {
            const size_t j = (size_t) (offset + n - 1);
            const double f = n * fFlex * std::sqrt (1.0 + str.inharmonicity * n * n);
            const double w = 2.0 * M_PI * f;
            const double sigma = internalDecay (f, str.fundamental, str.wound);
            const double sign = (n % 2 == 1) ? 1.0 : -1.0; // (-1)^(n+1)

            const cd pole = std::exp (cd (-sigma, w) * Ts);
            const cd poleOs = std::exp (cd (-sigma, w) * (Ts / os));
            const cd gain = pole * (Ts / w);          // impulse-invariant acceleration input
            const cd gainOs = poleOs * (Ts / os / w);

            m.omega[j] = (float) w;
            m.sigmaInternal[j] = (float) sigma;
            m.prFree[j] = (float) pole.real();
            m.piFree[j] = (float) pole.imag();
            m.prOs[j] = (float) poleOs.real();
            m.piOs[j] = (float) poleOs.imag();
            m.gOsr[j] = (float) gainOs.real();
            m.gOsi[j] = (float) gainOs.imag();

            // Bridge coupling. With the bridge end moving by u_b, mode n feels an
            // acceleration beta * u_b'' (beta = -2 (-1)^(n+1) / (n pi)), and the string pushes
            // on the bridge with T * dy/dx = T * (-1)^(n+1) * n pi / L per unit of mode amplitude.
            const double beta = -2.0 * sign / (n * M_PI);
            m.beta[j] = (float) beta;
            m.bgr[j] = (float) (gain.real() * beta);
            m.bgi[j] = (float) (gain.imag() * beta);
            m.fk[j] = (float) (str.tension * sign * n * M_PI / L);

            // Hammer: mode shape at the striking point, averaged over the felt's width
            const double x = n * M_PI * key.hammerWidth / (2.0 * L);
            const double width = std::abs (x) < 1e-9 ? 1.0 : std::sin (x) / x;
            m.phi[j] = (float) (std::sin (n * M_PI * strike) * width);
        }
        offset += plan.count;
    }

    setDamperStrength (damperStrength);
}

void ModalPiano::setDamperStrength (double scale)
{
    if (! m.omega.empty() && std::abs (scale - damperStrength) < 1e-3 && ! m.prDamp.empty() && m.prDamp[0] != 0.0f)
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
            for (int j = s.first; j < s.first + s.count; ++j)
            {
                m.pr[(size_t) j] = m.prDamp[(size_t) j];
                m.pi[(size_t) j] = m.piDamp[(size_t) j];
            }
}

void ModalPiano::buildSoundboard()
{
    // A grand soundboard has a roughly constant modal density (about one mode every 30 Hz).
    // Each mode is a damped resonator driven by the total bridge force. Mode masses are set
    // so that, averaged over frequency, the modes add about 1.5e-4 m/s/N to the bridge's
    // admittance (with the resistive part, about 3e-4 in all), in line with measured bridge
    // admittances (Wogram 1980; Conklin 1996; Giordano 1998). The low modes are sparse and
    // stand out as resonances, up to 1.5e-3 m/s/N.
    Rng rng (1234);
    std::vector<double> freqs;
    for (double f = 48.0; f < std::min (7000.0, 0.45 * sampleRate);)
    {
        freqs.push_back (f * (1.0 + 0.12 * (rng.next() - 0.5)));
        f += 33.0 * (0.6 + 0.8 * rng.next());
    }

    const size_t n = freqs.size();
    for (auto* v : { &board.wr, &board.wi, &board.qr, &board.qi, &board.hr, &board.hi, &board.omega, &board.sigma,
                     &board.radLeft, &board.radRight })
        v->assign (n, 0.0f);

    const double lossFactor = 0.02;
    const double meanBoardAdmittance = 1.5e-4;
    for (size_t i = 0; i < n; ++i)
    {
        const double f = freqs[i];
        const double w = 2.0 * M_PI * f;
        const double sigma = M_PI * f * lossFactor;
        const double meanSpacing = 2.0 * M_PI * 33.0;   // rad/s between modes
        const double peakAdmittance = std::min (1.5e-3, meanBoardAdmittance * meanSpacing / (M_PI * sigma));
        const double modalMass = 1.0 / (2.0 * sigma * peakAdmittance);

        const cd pole = std::exp (cd (-sigma, w) * Ts);
        const cd gain = pole * (Ts / (w * modalMass));
        board.omega[i] = (float) w;
        board.sigma[i] = (float) sigma;
        board.qr[i] = (float) pole.real();
        board.qi[i] = (float) pole.imag();
        board.hr[i] = (float) gain.real();
        board.hi[i] = (float) gain.imag();

        // Each mode radiates with its own strength and polarity towards each ear
        const double l = (0.4 + 0.6 * rng.next()) * (rng.next() < 0.5 ? -1 : 1);
        const double r = (0.4 + 0.6 * rng.next()) * (rng.next() < 0.5 ? -1 : 1);
        board.radLeft[i] = (float) l;
        board.radRight[i] = (float) r;
    }

    // Overall level: a mezzo-forte middle C peaks around -20 dBFS, fortissimo around -10
    outputGain = 45.0;
}

void ModalPiano::reset()
{
    std::fill (m.zr.begin(), m.zr.end(), 0.0f);
    std::fill (m.zi.begin(), m.zi.end(), 0.0f);
    std::fill (board.wr.begin(), board.wr.end(), 0.0f);
    std::fill (board.wi.begin(), board.wi.end(), 0.0f);
    lastBridgeVelocity = bridgeAcceleration = 0;
    resistiveLpLeft = resistiveLpRight = 0;
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
        s.damped = ! keys[(size_t) s.key].hasDamper ? false : true;
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
    for (int j = s.first + keep; j < s.first + s.count; ++j)
        m.zr[(size_t) j] = m.zi[(size_t) j] = 0.0f;
    s.activity = a;
}

void ModalPiano::updateDamping (int keyIndex)
{
    const auto& key = keys[(size_t) keyIndex];
    auto& ks = keyState[(size_t) keyIndex];
    const bool damped = key.hasDamper && ! ks.down && ! sustain;

    for (int i = ks.firstString; i < ks.firstString + ks.numStrings; ++i)
    {
        auto& s = strings[(size_t) i];
        if (s.damped != damped || s.activity == Activity::Off)
        {
            s.damped = damped;
            const auto& prSrc = damped ? m.prDamp : m.prFree;
            const auto& piSrc = damped ? m.piDamp : m.piFree;
            std::copy (prSrc.begin() + s.first, prSrc.begin() + s.first + s.count, m.pr.begin() + s.first);
            std::copy (piSrc.begin() + s.first, piSrc.begin() + s.first + s.count, m.pi.begin() + s.first);
        }
        // Undamped strings listen to the bridge
        if (! damped && s.activity == Activity::Off)
            setActivity (s, Activity::Sympathetic);
    }
}

void ModalPiano::setKeyDown (int midiNote, bool down)
{
    const int k = midiNote - kLowestKey;
    if (k < 0 || k >= (int) keys.size())
        return;
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

    // The hammer is never perfectly level, so it meets the strings of a unison a few
    // hundredths of a millimetre apart. This is what excites the slow "aftersound".
    Rng rng ((uint32_t) midiNote * 977u + 13u);
    double nearest = 1e9;
    for (int i = 0; i < ks.numStrings; ++i)
    {
        auto& s = strings[(size_t) (ks.firstString + i)];
        setActivity (s, Activity::Full);
        ks.hammer.gap[(size_t) i] = 0.04e-3 * rng.next();
        double u = 0;
        for (int j = s.first; j < s.first + s.count; ++j)
            u += m.phi[(size_t) j] * m.zi[(size_t) j];
        nearest = std::min (nearest, u + ks.hammer.gap[(size_t) i]);
    }

    ks.hammer.active = true;
    ks.hammer.position = nearest; // touching the first string
    ks.hammer.velocity = std::max (0.05, hammerVelocity);
    ks.hammer.samplesSinceStrike = 0;
    ks.hammer.samplesApart = 0;
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
    return m.omega[(size_t) (s.first + k - 1)] / (2.0 * M_PI);
}

void ModalPiano::checkActivity()
{
    // Retire strings that have gone quiet, and rebuild the work lists
    constexpr float silent = 1e-18f; // (1 nm)^2 displacement
    fullList.clear();
    hammerList.clear();
    activeStrings = 0;

    for (size_t k = 0; k < keys.size(); ++k)
    {
        auto& ks = keyState[k];
        if (ks.hammer.active)
            hammerList.push_back ((int) k);

        for (int i = ks.firstString; i < ks.firstString + ks.numStrings; ++i)
        {
            auto& s = strings[(size_t) i];
            if (s.activity == Activity::Off)
                continue;

            if (! ks.hammer.active)
            {
                const int splitAt = s.first + (s.activity == Activity::Full ? s.sympatheticCount : 0);
                float low = 0, high = 0;
                for (int j = s.first; j < splitAt; ++j)
                    low = std::max (low, m.zr[(size_t) j] * m.zr[(size_t) j] + m.zi[(size_t) j] * m.zi[(size_t) j]);
                const int end = s.first + (s.activity == Activity::Full ? s.count : s.sympatheticCount);
                for (int j = splitAt; j < end; ++j)
                    high = std::max (high, m.zr[(size_t) j] * m.zr[(size_t) j] + m.zi[(size_t) j] * m.zi[(size_t) j]);

                if (s.damped && std::max (low, high) < silent)
                    setActivity (s, Activity::Off);
                else if (! s.damped && s.activity == Activity::Full && high < silent && low < silent * 1e4f)
                    setActivity (s, Activity::Sympathetic);
            }

            if (s.activity == Activity::Full)
                ++activeStrings;
            if (s.activity != Activity::Off && ! ks.hammer.active)
                fullList.push_back (i);
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

    bool touching = false;
    for (int step = 0; step < os; ++step)
    {
        double hammerForce = 0;
        for (int i = 0; i < ks.numStrings; ++i)
        {
            auto& s = strings[(size_t) (ks.firstString + i)];
            const size_t a = (size_t) s.first, b = a + (size_t) s.count;

            double u = 0;
            for (size_t j = a; j < b; ++j)
                u += m.phi[j] * m.zi[j];

            const double compression = h.position - u - h.gap[(size_t) i];
            const double f = compression > 0 ? q0 * std::pow (compression * 1000.0, p) : 0.0;
            hammerForce += f;
            touching = touching || f > 0;

            const float hammerAcc = (float) (f / s.mass);
            const float bridgeAcc = (float) ub;
            for (size_t j = a; j < b; ++j)
            {
                const float acc = m.beta[j] * bridgeAcc + m.phi[j] * hammerAcc;
                const float zr = m.zr[j], zi = m.zi[j];
                m.zr[j] = m.prOs[j] * zr - m.piOs[j] * zi + m.gOsr[j] * acc;
                m.zi[j] = m.prOs[j] * zi + m.piOs[j] * zr + m.gOsi[j] * acc;
            }
        }

        // Hammer flight (symplectic Euler)
        h.velocity -= hammerForce / key.hammerMass * dt;
        h.position += h.velocity * dt;
    }

    // Force on the bridge from this key's strings
    for (int i = 0; i < ks.numStrings; ++i)
    {
        const auto& s = strings[(size_t) (ks.firstString + i)];
        double f = 0;
        for (int j = s.first; j < s.first + s.count; ++j)
            f += m.fk[(size_t) j] * m.zi[(size_t) j];
        forceSum += f;
        forceLeft += f * s.panLeft;
        forceRight += f * s.panRight;
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
        const float a = (float) bridgeAcceleration;
        double forceSum = 0, forceLeft = 0, forceRight = 0;

        // Strings without a hammer on them: one step per sample
        for (int index : fullList)
        {
            const auto& s = strings[(size_t) index];
            const int count = s.activity == Activity::Full ? s.count : s.sympatheticCount;
            const int j = s.first;
            const float f = stepModes (zr + j, zi + j, pr + j, pi + j, bgr + j, bgi + j, fk + j, count, a);
            forceSum += f;
            forceLeft += f * s.panLeft;
            forceRight += f * s.panRight;
        }

        // Strings under a hammer: oversampled
        for (size_t i = 0; i < hammerList.size(); ++i)
            runHammerKey (hammerList[i], bridgeAcceleration, forceSum, forceLeft, forceRight);

        // Soundboard
        double boardVelocity = 0, radiatedLeft = 0, radiatedRight = 0;
        stepBoard (board.wr.data(), board.wi.data(), board.qr.data(), board.qi.data(), board.hr.data(),
                   board.hi.data(), board.omega.data(), board.sigma.data(), board.radLeft.data(),
                   board.radRight.data(), (int) boardModes, (float) forceSum, boardVelocity, radiatedLeft,
                   radiatedRight);

        // Bridge motion, fed back to every string on the next sample
        const double bridgeVelocity = bridgeResistance * forceSum + boardVelocity;
        bridgeAcceleration = (bridgeVelocity - lastBridgeVelocity) * sampleRate;
        lastBridgeVelocity = bridgeVelocity;

        // Sound: soundboard modes plus the direct (dense high-frequency) bridge path
        resistiveLpLeft = forceLeft * bridgeResistance + lpCoeff * (resistiveLpLeft - forceLeft * bridgeResistance);
        resistiveLpRight = forceRight * bridgeResistance + lpCoeff * (resistiveLpRight - forceRight * bridgeResistance);
        const double outL = radiatedLeft + resistiveLpLeft;
        const double outR = radiatedRight + resistiveLpRight;

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
