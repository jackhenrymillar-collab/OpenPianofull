/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace openpiano::piano
{

// Physical parameters of one string.
struct StringParams
{
    double f0;  // Fundamental frequency [Hz]
    double L;   // Length [m]
    double rho; // Linear density [kg/m]
    double S;   // Cross-sectional area [m^2]
    double E;   // Young's modulus [N/m^2]
    double b1;  // First damping coefficient
    double b2;  // Second damping coefficient
};

// Physical parameters of one hammer.
struct HammerParams
{
    double Mh;      // Mass [kg]
    double p;       // Felt stiffness exponent
    double bH;      // Felt damping
    double K;       // Felt stiffness [N/m^p]
    double a;       // Striking point, as a fraction of the string length
    double gMeters; // Width of the contact area [m]
};

// A lossy stiff string struck by a felt hammer, solved with finite differences
// (Chaigne & Askenfelt 1994; Saitis 2008).
//
// This is a straight port of OpenPiano's original PianoString/Hammer code. The arithmetic is
// kept exactly as upstream (including its single-precision powf calls) so the regression
// test can prove that it renders bit-for-bit the same audio. Model improvements belong in
// later changes, made on purpose and measured against that baseline.
class StringModel
{
public:
    StringModel (int sampleRate, const StringParams& string, const HammerParams& hammer);

    // Strikes the string with the hammer moving at velocity [m/s].
    void hit (double hammerVelocity);

    // Raises the damping (key released, no pedal). strength 1 is the upstream damper.
    void damp (double strength = 1.0);

    // Restores the string's own damping (damper lifted).
    void undamp();

    // Scales the hammer felt stiffness. 1 is the nominal value.
    void setHammerHardness (double scale) { K = hammer.K * scale; }

    double nextSample();

    void reset();

    bool isActive() const { return active; }
    double fundamental() const { return string.f0; }
    uint32_t spatialSteps() const { return N; }

private:
    void computeCoefficients();
    void checkIfActive();
    double meanOfColumn (uint32_t column, uint32_t start, uint32_t stop, bool absolute) const;

    StringParams string;
    HammerParams hammer;

    int Fs;
    double Ts;

    // String damping in use (raised while damped)
    double b1, b2;

    // Derived string parameters
    double Ms, Te, c, rGyr, eps;
    uint32_t N = 0;          // Number of spatial steps
    uint32_t lenX = 0;
    double Xs = 0;           // Spatial step [m]

    // FD scheme coefficients (Chaigne's article)
    double D = 0, r = 0, a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0;

    // Hammer state and coefficients
    double K = 0;            // Felt stiffness in use
    double d1 = 0, d2 = 0, dF = 0;
    uint32_t XsContact = 0;
    std::vector<double> hammerMask;
    std::array<double, 4> eta {}; // Hammer displacement
    std::array<double, 4> Fh {};  // Hammer force

    // Displacement y[space][time], with a circular buffer of 4 time steps
    std::vector<std::array<double, 4>> y;
    uint8_t n0 = 3, n1 = 2, n2 = 1, n3 = 0;

    // Output pickup region
    uint32_t leftPickup = 0, rightPickup = 0;

    bool active = false;
    uint64_t activeCheckCounter = 0;
};

} // namespace openpiano::piano
