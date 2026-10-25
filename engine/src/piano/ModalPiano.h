/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include "ScaleDesign.h"

#include <algorithm>
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
// Each string vibrates in two directions: vertical (towards the soundboard, struck by the
// hammer) and horizontal (sideways, along the board). The horizontal motion pushes on the
// bridge much less, so it rings longer, and it leaks into the vertical through the bridge.
// Strings also stretch along their length (longitudinal waves): the stretch rises with the
// square of the vibration, so loud notes gain "phantom" partials and a metallic attack.
//
// Hammers: a mass on a nonlinear felt spring with memory (Stulov), F = Q0 (x^p - e * relaxed part), pressing
// on all of a key's strings. While a hammer touches the strings that key runs oversampled.
//
// Bridge and soundboard: every string pushes on the bridge, and the bridge's motion moves
// the end of every string. The bridge is a resistive part (the dense high-frequency
// soundboard response) plus a bank of soundboard and case modes. This one coupling gives:
//   - double decay: a key's 2 or 3 strings share the bridge, so their in-step motion drains
//     quickly into it while out-of-step motion lingers (Weinreich 1977), and the slow
//     horizontal motion adds a second, longer tail;
//   - sympathetic resonance: undamped strings (held keys, sustain pedal, top octaves with
//     no dampers, duplex segments behind the bridge) pick up energy from the bridge;
//   - the soundboard's own resonances colour and radiate the sound.
class ModalPiano
{
public:
    struct Settings
    {
        double unisonSpreadCents = 0.45;    // how far each outer unison string sits from the centre
        double hammerUnevenness = 0.04e-3;  // [m] how far apart a hammer meets a key's strings
        int oversampling = 8;               // during hammer contact
        bool horizontal = true;             // second (sideways) vibration direction
        bool longitudinal = true;           // lengthwise waves and phantom partials
        bool duplex = true;                 // undamped string segments behind the treble bridge
        double longitudinalCoupling = 0.12; // share of the lengthwise force that reaches the board (via bridge tilt)
        double crown = 1.0;                 // soundboard crown and downbearing (0 = flat, unloaded board)
        double feltMemory = 0.6;            // felt relaxation strength (0 = perfectly elastic)
        double feltRelaxTime = 0.6e-3;      // [s] how quickly squeezed felt relaxes
        double actionNoise = 1.0;           // key, hammer and damper noises (0 = silent action)
        double outerFelt = 0.03;            // stiffness of the soft outer felt, relative to the core
        double outerFeltScale = 1.0;        // scales how deep the outer felt is
        double farEnd = 1.0;                // agraffe / pressure-bar end moves with the iron frame
        double nonlinear = 1.0;             // instantaneous tension feedback on every partial
        double tensionGlide = 0.0;          // loud notes start sharp as the tension spike stretches the string
        double tensionLoss = 0.025;         // and lose their high partials faster while it lasts
        double air = 1.0;                   // air trapped under the board, between the rim and the floor
        double lid = 1.0;                   // lid opening (0 = closed, 1 = fully open)
        double age = 0.0;                   // 0 = new piano, 1 = old and worn (see prepare)
        double shankFlex = 1.0;             // hammer shank bending and head rocking (0 = rigid hammer)
        double softLift = 0.35;             // lifts the volume of soft notes without changing their tone (0 = off)
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
    void setActionNoise (double level) { settings.actionNoise = level; }

    // Adds nothing: overwrites left and right
    void process (float* left, float* right, int numSamples);

    int activeStringCount() const { return activeStrings; }
    int lowestNote() const { return kLowestKey; }
    int highestNote() const { return kHighestKey; }
    const std::vector<KeySpec>& scale() const { return keys; }
    const Settings& currentSettings() const { return settings; }

    // Frequency of a string's k-th partial as simulated (for tests)
    double partialFrequency (int midiNote, int stringIndex, int k) const;

private:
    enum class Activity : uint8_t { Off, Sympathetic, Full };
    enum class Kind : uint8_t { Vertical, Horizontal, Duplex, Front };

    struct String
    {
        int key;              // index into keys
        Kind kind = Kind::Vertical;
        int count = 0;            // modes in all
        int sympatheticCount = 0; // modes kept while only resonating (the lowest ones)
        int first = 0;            // where the lowest sympatheticCount modes sit in the mode arrays
        int highFirst = 0;        // where the rest sit
        float mass = 0;           // modal mass mu L / 2
        float panLeft = 0, panRight = 0;
        float hammerShare = 0;    // horizontal strings: fraction of the hammer force they feel
        int horizontal = -1;      // vertical strings: index of the same string's horizontal motion
        int longFirst = 0, longCount = 0, longPartials = 0; // longitudinal modes
        float tensionCoef = 0;    // EA/4 (pi/L)^2, for the average tension change
        mutable float tensionRise = 0; // latest tension rise / rest tension
        float restTension = 1;
        float outScale = 1;       // loudness lift for soft notes (bridge force up, bridge input down)
        int zone = 0;             // bridge zone the string presses on
        bool damped = true;
        Activity activity = Activity::Off;
    };

    struct Hammer
    {
        bool active = false;
        double position = 0, velocity = 0; // felt tip, moving towards the strings is positive
        double headPos = 0, headVel = 0;   // the head's wooden core, which rocks against the tip
        double shankPos = 0, shankVel = 0; // the shank, which bends under the head
        int samplesSinceStrike = 0, samplesApart = 0;
        std::array<double, 3> gap {};          // hammer not perfectly level across the unison
        std::array<double, 3> lastCompression {};
        std::array<double, 3> feltHistory {};  // relaxing part of the felt force, per string
        std::array<double, 3> forceSum {};     // per string, over the current sample
    };

    struct KeyState
    {
        int firstString = 0, numStrings = 0; // vertical strings
        int firstHorizontal = -1;            // their horizontal motion, same order
        int duplexString = -1;
        bool tensionApplied = false;
        long long releasedAt = -1;   // sample clock when the key was let go
        int damperPending = 0;       // samples until the rising key lets the damper land
        bool down = false;
        Hammer hammer;
    };

    void buildModes();
    int addString (const StringSpec& spec, int key, Kind kind, double strikeRatio, double hammerWidth,
                   int maxModes, double maxFrequency, uint32_t seed);
    void addLongitudinal (String& s, const StringSpec& spec);
    void buildSoundboard();
    void updateDamping (int keyIndex);
    void setActivity (String& s, Activity a);
    void runHammerKey (int keyIndex, double bridgeAcceleration, double& forceSum, double& forceLeft, double& forceRight);
    double runLongitudinal (const String& s);
    void checkActivity();
    void applyTension();
    int tensionCounter = 0;
    long long sampleClock = 0;

    // Each string's modes are stored in two parts: the low ones it keeps while only
    // resonating sit in one region, back to back with every other string's, and the rest
    // in another. So resonating strings form long runs that can be stepped in one go.
    static int modeIndex (const String& s, int k)
    {
        return k < s.sympatheticCount ? s.first + k : s.highFirst + k - s.sympatheticCount;
    }
    template <typename Fn> static void forModes (const String& s, int keep, Fn&& fn)
    {
        const int low = std::min (keep, s.sympatheticCount);
        for (int j = s.first; j < s.first + low; ++j)
            fn ((size_t) j);
        for (int j = s.highFirst; j < s.highFirst + keep - low; ++j)
            fn ((size_t) j);
    }
    void relayoutModes();

    template <typename Fn> void forEachStringConst (const KeyState& ks, Fn&& fn) const
    {
        for (int i = 0; i < ks.numStrings; ++i)
        {
            fn (strings[(size_t) (ks.firstString + i)]);
            if (ks.firstHorizontal >= 0)
                fn (strings[(size_t) (ks.firstHorizontal + i)]);
        }
    }
    template <typename Fn> void forEachString (const KeyState& ks, Fn&& fn)
    {
        for (int i = 0; i < ks.numStrings; ++i)
        {
            fn (strings[(size_t) (ks.firstString + i)]);
            if (ks.firstHorizontal >= 0)
                fn (strings[(size_t) (ks.firstHorizontal + i)]);
        }
    }

    double sampleRate = 48000, Ts = 1.0 / 48000;
    Settings settings;
    double hardness = 1.0;
    double trebleFirmness = 2.2;
    double ageFactor = 0;
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
        std::vector<float> gr, gi;          // acceleration input gain, one sample
        std::vector<float> bgr, bgi;        // bridge-acceleration input gain (gain * beta)
        std::vector<float> fk;              // contribution to the bridge force
        std::vector<float> alt;             // (-1)^(n+1): turns bridge-end coupling into far-end coupling
        std::vector<float> beta;            // bridge coupling weight
        std::vector<float> phi;             // mode shape at the hammer (incl. hammer width)
        std::vector<float> prOs, piOs, gOsr, gOsi; // oversampled pole and input gain
        std::vector<float> sigmaInternal, omega;

        void clear();
        size_t size() const { return zr.size(); }
    } m;

    // Longitudinal modes
    struct Longitudinal
    {
        std::vector<float> zr, zi, pr, pi, gr, gi;
        std::vector<float> drive;   // acceleration per unit of the pair sum
        std::vector<float> out;     // bridge force per unit of displacement
        std::vector<float> invW2;   // 1 / omega^2, for the slowly varying part
        std::vector<int> order;     // m
    } lm;
    std::vector<float> scratchU;

    // Soundboard and case modes
    struct Board
    {
        std::vector<float> wr, wi, qr, qi, hr, hi, omega, sigma, radLeft, radRight;
        std::vector<float> phi;    // mode shape at each bridge zone, zone by zone
        std::vector<float> scratchF, scratchV;
        std::vector<float> volume; // net volume each mode pushes (drives the air under the board)
    } board;

    // Air under the board: a few broad cavity resonances driven by the board's volume motion.
    // They radiate (mostly the same into both ears) and press back on the board.
    struct AirMode { float b0, a1, a2, x1, x2, y1, y2, gain; };
    std::vector<AirMode> airModes;
    double airPressure = 0;
    // Lid: an early, treble-heavy reflection towards the listener, a little later on one side
    std::vector<float> lidBufLeft, lidBufRight;
    int lidPos = 0, lidDelayLeft = 0, lidDelayRight = 0;
    double lidHpLeft = 0, lidHpRight = 0, lidHpInLeft = 0, lidHpInRight = 0, lidHpCoeff = 0, lidLpLeft = 0, lidLpRight = 0, lidLpCoeff = 0;
    void buildBody();

    // Action noises: short knocks from the keys, hammers and dampers, ringing a small bank
    // of wooden resonances (key, keybed, hammer shank) and reaching the board through the frame
    struct ActionHit { int delay; int type; float amount, panLeft, panRight; };
    std::vector<ActionHit> actionHits;
    struct ActionBank
    {
        std::array<float, 8> zr {}, zi {}, pr {}, pi {};
        std::array<std::array<float, 8>, 3> weight {}; // per hit type
        float panLeft = 0.7f, panRight = 0.7f;
    } action;
    void buildAction();
    void queueAction (int type, int delaySamples, double amount, int keyIndex);
    void runAction (double& force, double& left, double& right);

    double bridgeResistance = 0.9e-4;      // m/s per N: the dense high-frequency part of the soundboard
    double horizontalResistance = 0.6e-4;  // the bridge moves less easily sideways
    double crossResistance = 0.36e-4;      // vertical force -> sideways motion, and back
    // The bridge is split into zones along its length; each key's strings press on their own
    // zone, and every board mode moves each zone by its own amount
    static constexpr int kZones = 8;
    std::vector<int> keyZone;
    std::array<double, kZones> zoneForce {}, zoneAcceleration {}, lastZoneVelocity {}, zoneBoardVelocity {};
    double lastBridgeVelocity = 0, bridgeAcceleration = 0;
    double lastHorizontalVelocity = 0, horizontalAcceleration = 0;
    double resistiveLpLeft = 0, resistiveLpRight = 0, lpCoeff = 0;
    double dcLeft = 0, dcRight = 0, dcInLeft = 0, dcInRight = 0, dcCoeff = 0;
    double outputGain = 1.0;
    double airLoad = 5.0; // N per unit of air pressure on the board

    int activityCounter = 0;
    int activeStrings = 0;
    std::vector<int> fullList, hammerList, longList, alwaysOn;
    std::vector<int> frameList; // strings whose far end is worked by the frame this block

    // The iron frame at the strings' far ends (agraffes, pressure bar): a stiff, heavy,
    // barely damped structure. Its few modes are pushed by every sounding string's far end
    // and move those ends back, coupling the strings a second way.
    struct Frame { std::vector<float> wr, wi, qr, qi, hr, hi, omega, sigma, rad; double accel = 0, lastVel = 0; } frame;
    void buildFrame();
    void runFrame (double& left, double& right);

    // Neighbouring strings whose modes sit back to back are stepped in one go
    struct Run { int first, count; float panLeft, panRight; bool side; int zone; };
    std::vector<Run> runs;
    struct Segment { int first, count; float panLeft, panRight; bool side; int zone; };
    std::vector<Segment> segmentScratch;
};

} // namespace openpiano::piano
