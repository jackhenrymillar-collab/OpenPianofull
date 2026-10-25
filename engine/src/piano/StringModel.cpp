/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#include "StringModel.h"

#include <algorithm>
#include <cmath>

namespace openpiano::piano
{

// NOTE: the powf() calls below are deliberate. Upstream computed these terms in single
// precision, and we keep that until the string model is reworked so renders stay identical.

StringModel::StringModel (int sampleRate, const StringParams& s, const HammerParams& h)
    : string (s), hammer (h), Fs (sampleRate), Ts (1. / sampleRate), b1 (s.b1), b2 (s.b2)
{
    // Hammer felt coefficients (Saitis, Eq. 4.21)
    K = hammer.K;
    d1 = 2 / (1 + hammer.bH * Ts / (2 * hammer.Mh));
    d2 = (-1 + hammer.bH * Ts / (2 * hammer.Mh)) / (1 + hammer.bH * Ts / (2 * hammer.Mh));
    dF = (-powf (Ts, 2) / hammer.Mh) / (1 + hammer.bH * Ts / (2 * hammer.Mh));

    // Derived physical parameters
    Ms = string.rho * string.L;
    Te = string.rho * powf (string.L, 2) * 4 * powf (string.f0, 2);
    c = std::sqrt (Te / string.rho);
    rGyr = string.S / 2;
    eps = powf (rGyr, 2) * ((string.E * string.S) / (Te * powf (string.L, 2)));

    // Largest number of spatial steps that avoids aliasing (Chaigne, Eq. 11-12)
    const double gamma = Fs / (2 * string.f0);
    N = (uint32_t) std::floor (std::sqrt ((-1 + std::sqrt (1 + 16 * eps * powf (gamma, 2))) / (8 * eps)));
    lenX = N;
    Xs = string.L / N;

    // Hammer contact window: a Hann window centred on the striking point
    const double xContact = hammer.a * string.L;
    XsContact = (uint32_t) std::round (xContact / Xs);
    const double g = std::ceil (hammer.gMeters * N / string.L);
    const auto gLen = (uint32_t) g;
    hammerMask.assign (lenX, 0.0);
    const int start = (int) floorf (XsContact - (g / 2)) + 1;
    for (uint32_t k = 0; k < gLen; ++k)
        hammerMask[(size_t) start + k] = 0.5 * (1 - std::cos (2 * M_PI * k / (gLen - 1)));

    y.assign (lenX + 2, { 0.0, 0.0, 0.0, 0.0 });

    // Pickup: the mean over a few points mirrored from the striking point
    const uint32_t nSpace = std::min ((uint32_t) 13, ((N - 1) | 0x1));
    const uint32_t xsSound = N - XsContact;
    leftPickup = xsSound - (nSpace - 1) / 2;
    rightPickup = xsSound + (nSpace - 1) / 2;

    computeCoefficients();
}

void StringModel::computeCoefficients()
{
    // FD scheme coefficients (Chaigne, Eq. 10)
    r = c * Ts / Xs;
    const double rSqr = r * r;
    const double NSqr = (double) N * N;
    D = 1 + b1 * Xs + 2 * b2 / Ts;
    a1 = (2 - 2 * rSqr + b2 / Ts - 6 * eps * NSqr * rSqr) / D;
    a2 = (-1 + b1 * Ts + 2 * b2 / Ts) / D;
    a3 = (rSqr * (1 + 4 * eps * NSqr)) / D;
    a4 = (b2 / Ts - eps * NSqr * rSqr) / D;
    a5 = (-b2 / Ts) / D;
}

void StringModel::reset()
{
    for (auto& column : y)
        column = { 0.0, 0.0, 0.0, 0.0 };
    eta = {};
    Fh = {};
    n0 = 3; n1 = 2; n2 = 1; n3 = 0;
    active = false;
    activeCheckCounter = 0;
    undamp();
}

double StringModel::meanOfColumn (uint32_t column, uint32_t start, uint32_t stop, bool absolute) const
{
    double sum = 0;
    for (uint32_t i = start; i < stop; ++i)
        sum += absolute ? std::fabs (y[i][column]) : y[i][column];
    return sum / (stop - start);
}

void StringModel::checkIfActive()
{
    // Sum of the mean absolute displacement over the 4 buffered time steps.
    // The 1 micrometre threshold is empirical.
    double displacement = 0.0;
    displacement += meanOfColumn (n0, 2, lenX - 3, true);
    displacement += meanOfColumn (n1, 2, lenX - 3, true);
    displacement += meanOfColumn (n2, 2, lenX - 3, true);
    displacement += meanOfColumn (n3, 2, lenX - 3, true);
    active = displacement > 0.000001;
}

void StringModel::hit (double hammerVelocity)
{
    active = true;
    activeCheckCounter = 0;
    undamp();

    // Step the time indices back by one, so the next nextSample() call lands on n0 again.
    n3 = (n3 - 1) & 0x3;
    n2 = (n2 - 1) & 0x3;
    n1 = (n1 - 1) & 0x3;
    n0 = (n0 - 1) & 0x3;

    // The hammer starts in contact with the string, moving at hammerVelocity
    eta[n3] = 0;
    eta[n2] = 0;
    eta[n1] = 0;
    eta[n0] = hammerVelocity * Ts;

    if (eta[n0] < y[XsContact][n0]) // Chaigne, Eq. 21
        Fh[n0] = 0.0f;
    else
        Fh[n0] = K * powf (eta[n0] - y[XsContact][n0], hammer.p); // Chaigne, Eq. 20
}

void StringModel::undamp()
{
    b1 = string.b1;
    b2 = string.b2;
    computeCoefficients();
}

void StringModel::damp (double strength)
{
    b1 = 0.2 * strength;
    b2 = 6.25e-6 * strength;
    computeCoefficients();
}

double StringModel::nextSample()
{
    // Skip strings whose displacement is negligible. Checked every 16384 samples.
    if (++activeCheckCounter > 0x4000)
    {
        activeCheckCounter = 0;
        checkIfActive();
    }
    if (! active)
        return 0;

    n3 = (n3 + 1) & 0x3;
    n2 = (n2 + 1) & 0x3;
    n1 = (n1 + 1) & 0x3;
    n0 = (n0 + 1) & 0x3;

    // String displacement (Chaigne, Eq. 10)
    const double forceScale = Ts * Ts * N * Fh[n1];
    for (uint32_t i = 2; i < lenX - 3; ++i)
    {
        y[i][n0] = a1 * y[i][n1] + a2 * y[i][n2]
                 + a3 * (y[i + 1][n1] + y[i - 1][n1])
                 + a4 * (y[i + 2][n1] + y[i - 2][n1])
                 + a5 * (y[i + 1][n2] + y[i - 1][n2] + y[i][n3])
                 + (forceScale * hammerMask[i]) / Ms;
    }

    // Simplified boundaries with perfect reflection (Chaigne, Eq. 23)
    const uint32_t end = lenX + 1;
    y[0][n0] = -y[2][n0];
    y[end][n0] = -y[end - 2][n0];

    // Hammer displacement including felt damping (Saitis, Eq. 4.21)
    eta[n0] = d1 * eta[n1] + d2 * eta[n2] + dF * Fh[n1];

    // Hammer force, zero when the hammer has left the string (Chaigne, Eq. 20-21)
    if (eta[n0] < y[XsContact][n0])
        Fh[n0] = 0.0f;
    else
        Fh[n0] = K * powf (eta[n0] - y[XsContact][n0], hammer.p);

    return meanOfColumn (n0, leftPickup, rightPickup, false);
}

} // namespace openpiano::piano
