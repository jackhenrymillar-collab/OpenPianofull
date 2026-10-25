/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include "ScaleDesign.h"

#include <array>
#include <cstdint>
#include <vector>

namespace openpiano::piano
{

// The whole piano as one coupled physical system:
//
//   hammers --(nonlinear felt)--> strings --(bridge force)--> bridge / soundboard --> sound
//                                    ^                              |
//                                    +------(bridge motion)---------+
//
// Strings: each string solves the lossy stiff-string equation that OpenPiano used, but in
// its vibration modes instead of on a spatial grid (modal synthesis). Mode k sits at
// k f sqrt(1 + B k^2) and decays at its own rate. There is no grid to run out of, so treble
// notes are as accurate as bass notes, at any sample rate.
//
// Hammers: a mass on a nonlinear felt spring, F = Q0 x^p, pressing on all of a key's
// strings. While a hammer touches the strings that key runs 8x oversampled.
//
// Bridge and soundboard: every string pushes on the bridge, and the bridge's motion moves
// the end of every string. The bridge is a resistive part (the dense high-frequency
// soundboard response) plus a bank of soundboard modes. This one coupling gives:
//   - double decay: a key's 2 or 3 strings share the bridge, so their in-step motion drains
//     quickly into it while out-of-step motion lingers (Weinreich 1977);
//   - sympathetic resonance: undamped strings (held keys, sustain pedal, top octaves with
//     no dampers) pick up energy from the bridge;
//   - the soundboard's own resonances colour and radiate the sound.
class ModalPiano
{
public:
    struct Settings
    {
        double unisonSpreadCents = 0.3; // how far apart unison strings are tuned
        int oversampling = 8;           // during hammer contact
    };

    ModalPiano();

    void prepare (double sampleRate, const Settings& settings);
    void reset();

    // velocity in m/s at the hammer
    void strike (int midiNote, double hammerVelocity);
    void setKeyDown (int midiNote, bool down);
    void setSustain (bool down);

    void setHammerHardness (double scale) { hardness = scale; }
    void setDamperStrength (double scale);

    // Adds nothing: overwrites left and right
    void process (float* left, float* right, int numSamples);

    int activeStringCount() const { return activeStrings; }
    int lowestNote() const { return kLowestKey; }
    int highestNote() const { return kHighestKey; }
    const std::vector<KeySpec>& scale() const { return keys; }

    // Frequency of a string's k-th partial as simulated (for tests)
    double partialFrequency (int midiNote, int stringIndex, int k) const;

private:
    enum class Activity : uint8_t { Off, Sympathetic, Full };

    struct String
    {
        int key;              // index into keys
        int first, count;     // range in the mode arrays
        int sympatheticCount; // modes kept while only resonating
        float mass;           // modal mass mu L / 2
        float panLeft, panRight;
        bool damped = true;
        Activity activity = Activity::Off;
    };

    struct Hammer
    {
        bool active = false;
        double position = 0, velocity = 0; // moving towards the strings is positive
        int samplesSinceStrike = 0, samplesApart = 0;
        std::array<double, 3> gap {};      // hammer not perfectly level across the unison
    };

    struct KeyState
    {
        int firstString, numStrings;
        bool down = false;
        Hammer hammer;
    };

    void buildModes();
    void buildSoundboard();
    void updateDamping (int keyIndex);
    void setActivity (String& s, Activity a);
    void runHammerKey (int keyIndex, double bridgeAcceleration, double& forceSum, double& forceLeft, double& forceRight);
    void checkActivity();

    double sampleRate = 48000, Ts = 1.0 / 48000;
    Settings settings;
    double hardness = 1.0;
    double damperStrength = 1.0;
    bool sustain = false;

    std::vector<KeySpec> keys;
    std::vector<KeyState> keyState;
    std::vector<String> strings;

    // Per-mode data, all strings back to back (structure of arrays, float for SIMD)
    struct Modes
    {
        std::vector<float> zr, zi;          // complex state; zi is the displacement
        std::vector<float> pr, pi;          // per-sample pole in use
        std::vector<float> prFree, piFree;  // pole without damper
        std::vector<float> prDamp, piDamp;  // pole with damper
        std::vector<float> bgr, bgi;        // bridge-acceleration input gain
        std::vector<float> fk;              // contribution to the bridge force
        std::vector<float> beta;            // bridge coupling weight
        std::vector<float> phi;             // mode shape at the hammer (incl. hammer width)
        std::vector<float> prOs, piOs, gOsr, gOsi; // oversampled pole and input gain
        std::vector<float> sigmaInternal, omega;

        void resize (size_t n);
    } m;

    // Soundboard modes
    struct Board
    {
        std::vector<float> wr, wi, qr, qi, hr, hi, omega, sigma, radLeft, radRight;
    } board;

    double bridgeResistance = 1.5e-4;  // m/s per N: the dense high-frequency part of the soundboard
    double lastBridgeVelocity = 0, bridgeAcceleration = 0;
    double resistiveLpLeft = 0, resistiveLpRight = 0, lpCoeff = 0;
    double dcLeft = 0, dcRight = 0, dcInLeft = 0, dcInRight = 0, dcCoeff = 0;
    double outputGain = 1.0;

    int activityCounter = 0;
    int activeStrings = 0;
    std::vector<int> fullList, hammerList;
};

} // namespace openpiano::piano
