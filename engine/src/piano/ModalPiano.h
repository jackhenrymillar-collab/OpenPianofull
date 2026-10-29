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
#include <complex>
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
        double unisonSpreadCents = 1.0;      // how far each outer unison string sits from the centre
        double hammerUnevenness = 0.04e-3;  // [m] how far apart a hammer meets a key's strings
        int oversampling = 8;               // during hammer contact
        bool horizontal = true;             // second (sideways) vibration direction
        bool longitudinal = true;           // lengthwise waves and phantom partials
        bool duplex = true;                 // undamped string segments behind the treble bridge
        double longitudinalCoupling = 0.2;  // share of the lengthwise force that reaches the board (via bridge tilt); 0.2 fits the recorded bass's loud brightness
        int longitudinalOrders = 16;        // lengthwise modes per string
        int longitudinalPartials = 64;     // transverse partials that drive them
        double longitudinalLoss = 40.0;     // decay of free lengthwise ringing [1/s]
        double longitudinalResonantShare = 0.35; // of that, for the lengthwise resonances themselves
        double longitudinalLossSlope = 0.15; // and its rise with frequency [1/s per Hz]
        bool contactSpread = true;          // the felt's contact patch widens as it is squeezed
        double crown = 1.0;                 // soundboard crown and downbearing (0 = flat, unloaded board)
        double feltMemory = 0.95;           // felt relaxation strength (0 = perfectly elastic)
        double feltRelaxTime = 3.0e-3;     // [s] how quickly squeezed felt relaxes
        double actionNoise = 1.0;           // key, hammer and damper noises (0 = silent action)
        double outerFelt = 0.02;            // stiffness of the soft outer felt, relative to the core
        double outerFeltScale = 1.0;        // scales how deep the outer felt is
        double farEnd = 1.0;                // agraffe / pressure-bar end moves with the iron frame
        double nonlinear = 0.1;             // instantaneous tension feedback on every partial (stronger makes
                                            // loud mid notes flash bright and fade at once, a plucked start)
        double tensionGlide = 0.0;          // loud notes start sharp as the tension spike stretches the string
        double tensionLoss = 0.0;         // and lose their high partials faster while it lasts
        double air = 1.0;                   // air trapped under the board, between the rim and the floor
        double lid = 1.0;                   // lid opening (0 = closed, 1 = fully open)
        bool closeMics = false;             // listen from a close stereo mic pair (as on studio sample sets) instead of the player's seat
        double age = 0.0;                   // 0 = new piano, 1 = old and worn (see prepare)
        double shankFlex = 1.0;             // hammer shank bending and head rocking (0 = rigid hammer)
        double softLift = 0.0;              // lifts the volume of soft bass notes without changing their tone, fading to none at the top (0 = off)
    };

    ModalPiano();

    void prepare (double sampleRate, const Settings& settings);
    void reset();

    // velocity in m/s at the hammer
    void strike (int midiNote, double hammerVelocity);
    // releaseSpeed (0..1, 0.5 = an ordinary release) is how fast the key was let go: it sets
    // how soon the damper reaches the strings and how hard it lands
    void setKeyDown (int midiNote, bool down, double releaseSpeed = 0.5);
    void setSustain (bool down) { setSustainLevel (down ? 1.0 : 0.0); }
    // Sustain pedal position, 0 = up, 1 = fully down. The first part of the travel does
    // nothing, then the dampers lift gradually (half pedalling), free from about 0.6
    void setSustainLevel (double level);
    // Sostenuto (middle pedal): pressed, it holds up the dampers that are up at that moment
    // (keys down, or all of them with the sustain pedal down); later notes are damped as usual
    void setSostenuto (bool down);
    // Una corda (left pedal), 0 = up, 1 = fully down. Applies to notes struck from then on.
    void setSoftPedal (double amount) { softPedal = std::clamp (amount, 0.0, 1.0); }

    void setHammerHardness (double scale) { hardness = scale; }
    void setDamperStrength (double scale);
    void setActionNoise (double level) { settings.actionNoise = level; }
    void setRoom (double amount) { room = (float) std::clamp (amount, 0.0, 1.0); }

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
        float sympTopF = 0;       // frequency of the highest mode kept while only resonating [Hz]
        float strikeBase = 0;     // where the hammer meets it, as a fraction of the length
        float halfWidthOverL = 0; // hammer width / (2 L)
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
        double feltSet = 0;        // copied from the key at the strike
        std::array<double, 3> gap {};          // hammer not perfectly level across the unison
        std::array<double, 3> lastCompression {};
        std::array<double, 3> feltHistory {};  // relaxing part of the felt force, per string
        std::array<double, 3> forceSum {};     // per string, over the current sample
        std::array<double, 3> feltScale { 1, 1, 1 }; // each string meets slightly different felt
        double softFelt = 0;                   // una corda: struck with the softer, unworn felt
    };

    struct KeyState
    {
        int firstString = 0, numStrings = 0; // vertical strings
        int firstHorizontal = -1;            // their horizontal motion, same order
        int duplexString = -1;
        bool tensionApplied = false;
        long long releasedAt = -1;   // sample clock when the key was let go
        long long struckAt = -1;     // and when it was last struck
        double feltSet = 0;          // how far the felt is still squeezed from earlier blows (0-1)
        int damperPending = 0;       // samples until the rising key lets the damper land
        int damperRamp = 0;          // samples left while the landed felt settles into its full grip
        float releaseSpeed = 0.5f;   // how fast the key was let go, 0..1
        bool down = false;
        bool sostenuto = false;      // damper held up by the sostenuto pedal
        Hammer hammer;
    };

    void buildModes();
    int addString (const StringSpec& spec, int key, Kind kind, double strikeRatio, double hammerWidth,
                   int maxModes, double maxFrequency, uint32_t seed);
    void addLongitudinal (String& s, const StringSpec& spec);
    void buildSoundboard();
    void updateDamping (int keyIndex);
    // The pole a string's mode rests at now: free, damped, or part way while the felt settles
    std::complex<double> restingPole (const String& s, const KeyState& ks, size_t j) const;
    void setActivity (String& s, Activity a);
    void runHammerKey (int keyIndex, double bridgeAcceleration, double& forceSum, double& forceLeft, double& forceRight);
    double runLongitudinal (const String& s);
    void checkActivity();
    void applyTension();
    int tensionCounter = 0;
    long long sampleClock = 0;
    std::vector<double> tuneFactor;
    static constexpr std::array<double, 88> kTunerCents { 3.2, 2.2, 1.4, 0.6, 3.3, 1.2, 0.9, 1.2, 4.8, 4.1, 0.3, 2.7, 0.0, 1.0, 4.0, 2.7, 3.8, 1.5, 1.6, -0.2, -0.16, -3.83, 1.66, 1.94, -0.52, -1.19, -3.12, -3.79, -4.1, -1.3, 2.74, 2.81, 2.57, 0.67, -0.99, -1.97, -1.7, -0.77, -1.4, 1.6, -2.46, -1.1, 2.29, 2.2, 2.99, 2.83, 1.9, 2.06, 1.78, 1.5, -2.21, 1.32, -1.12, 1.02, -5.21, -5.17, -5.33, 1.64, 1.07, -7.17, -7.25, -8.36, 0.72, -0.8, -1.42, 0.77, -3.43, -1.05, 0.75, 0.39, -0.46, 0.21, 0.42, 0.25, -0.35, -0.11, 0.2, -0.93, -1.2, 1.31, -1.12, 0.46, -0.83, 0.17, 0.22, 0.79, -0.92, -0.33 };   // filled by the tuning pass   // the tuner's correction per key, for the body's pull on its pitch

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
    double trebleFirmness = 0.5;
    double ageFactor = 0;
    double damperStrength = 1.0;
    bool sustain = false;        // the pedal holds every damper clear of the strings
    double pedalLevel = 0;       // sustain pedal position, 0..1
    double pedalGrip = 1.0;      // how much of their grip the dampers keep (1 = pedal up)
    bool sostenutoDown = false;
    double softPedal = 0;

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
        std::vector<float> pt, wa;          // the same at a single point, and the width term (see runHammerKey)
        std::vector<float> prOs, piOs, gOsr, gOsi; // oversampled pole and input gain
        std::vector<float> sigmaInternal, omega;
        std::vector<float> damperShift;     // (damped frequency / free frequency) - 1: the damper felt and mass loading the string

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
        std::vector<float> wr, wi, qr, qi, hr, hi, direct, omega, sigma, radLeft, radRight;
        std::vector<float> phi;    // mode shape at each bridge zone, zone by zone
        std::vector<float> scratchF, scratchV;
        std::vector<float> phiH;
        std::vector<float> phiSum; // every zone's shape added up: a push spread evenly along the bridge   // sideways shapes at the zones (bridge tilt), same layout as phi
        std::vector<float> volume; // net volume each mode pushes (drives the air under the board)
    } board;

    // Air under the board: a few broad cavity resonances driven by the board's volume motion.
    // They radiate (mostly the same into both ears) and press back on the board.
    struct AirMode { float b0, a1, a2, x1, x2, y1, y2, gain; };
    std::vector<AirMode> airModes;
    double airPressure = 0;
    // Lid: an early, treble-heavy reflection towards the listener, a little later on one side
    void buildBody();
    // Above a few kHz the board is a dense tangle of overlapping modes: what the bridge
    // pushes comes out smeared over a millisecond or two and differently in each ear,
    // not as a click. Three allpass diffusers per ear model that.
    struct Diffuser { std::vector<float> buf; int pos = 0; float g = 0.7f; float process (float x)
        { const float d = buf[(size_t) pos]; const float y = -g * x + d; buf[(size_t) pos] = x + g * y; pos = (pos + 1) % (int) buf.size(); return y; } };
    std::array<Diffuser, 3> diffLeft, diffRight;
    // Simple biquad for the body's colouring filters
    struct Biquad { float b0 = 1, b1 = 0, b2 = 0, a1 = 0, a2 = 0, x1 = 0, x2 = 0, y1 = 0, y2 = 0;
        float process (float x) { const float y = b0 * x + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
                                  x2 = x1; x1 = x; y2 = y1; y1 = y; return y; } };
    std::vector<Biquad> highBoardLeft, highBoardRight; // board peaks and dips above the simulated modes
    // Early reflections of the room, as heard from the player's seat: floor, side walls,
    // ceiling and the wall behind, each a delayed, duller copy of the piano's sound
    struct Reflection { int delay; int fromRight; float gain, lpCoeff, state; bool toRight; };
    std::vector<Reflection> reflections;
    std::vector<float> roomBufLeft, roomBufRight;
    int roomPos = 0;
    float room = 0.35f;
    double beamLpLeft = 0, beamLpRight = 0, beamCoeff = 0; // high frequencies beam towards the open side

    // Action noises: short knocks from the keys, hammers and dampers, ringing a small bank
    // of wooden resonances (key, keybed, hammer shank) and reaching the board through the frame
    struct ActionHit { int delay; int type; float amount, panLeft, panRight; };
    std::vector<ActionHit> actionHits;
    struct ActionBank
    {
        static constexpr int N = 11;
        static constexpr int Types = 4;
        std::array<float, N> zr {}, zi {}, pr {}, pi {};
        std::array<std::array<float, N>, Types> weight {}; // per hit type
        float panLeft = 0.7f, panRight = 0.7f;
    } action;
    void buildAction();
    void queueAction (int type, int delaySamples, double amount, int keyIndex);
    void runAction (double& force, double& left, double& right);

    double dbgRad = 1, dbgDirect = 1; std::array<double, 8> dbgRZ { 1, 1, 1, 1, 1, 1, 1, 1 }; // TEMP
    double bridgeResistance = 0.9e-4;      // m/s per N: the dense high-frequency part of the soundboard
    double feltScatter = 0.12;            // how unevenly one hammer's felt meets the strings of its unison
    double frameKnockGain = 150.0;          // N per (m/s)^1.5: how hard the hammer blow shakes the iron plate
    // Side bearing: the strings zig-zag round the bridge pins, which lean outwards. Pressed
    // against a leaning pin, a string's sideways pull has an upward part, and the bridge's
    // up-and-down motion moves the string's end sideways by the same share (about tan of the
    // pins' lean). This couples the two directions of vibration at the bridge itself.
    double sideBearing = 0.15;
    // The bridge is split into zones along its length (three on the bass bridge, thirteen on the
    // treble bridge); each key's strings press on their own zone, and every board mode moves each
    // zone by its own amount
    static constexpr int kZones = 32, kBassZones = 4;   // about 2.4 keys per zone on the treble bridge, 5 on the bass bridge
    std::vector<int> keyZone;
    std::array<double, kZones> bridgeMassCoeff {}, bridgeMassState {};
    double actionSpread = 0;                    // the action's push, spread along the bridge
    std::array<uint8_t, kZones> zoneLive {};   // zones with a string being run this block
    std::array<uint8_t, kZones> zoneSideLive {};   // zones with a sideways string being run this block
    std::array<double, kZones> zoneSideForce {}, zoneSideAccel {}, lastZoneSideVel {};
    std::array<double, kZones> zoneSideBoardVel {};   // the board's own sideways motion at each zone
    // How many of the board's modes (sorted by frequency) each zone needs this block: all of them
    // under a struck string, but only those up to the highest partial kept by strings that are
    // merely resonating. A resting string's few low partials cannot push or hear the modes above.
    std::array<int, kZones> zoneModesV {}, zoneModesS {};
    std::array<double, kZones> highAdmittance {}, highPassState {}, highPassState2 {};   // the many-mode board above the solved modes
    double highPassCoeff = 0, highRadiation = 0, highLeftState = 0, highRightState = 0, highLeftState2 = 0, highRightState2 = 0;
    // Where a zone sits along the bridges on the scale of the earlier eight-zone layout (0..7):
    // 0..1 the bass bridge, 2..7 the treble bridge from its low end to the top
    static double layoutPosition (int zone);
    std::array<double, kZones> zoneCoupling {};
    std::array<double, kZones> zoneForce {}, zoneAcceleration {}, lastZoneVelocity {}, zoneBoardVelocity {};
    double lastBridgeVelocity = 0, bridgeAcceleration = 0;
    double lastHorizontalVelocity = 0, horizontalAcceleration = 0;
    double resistiveLpLeft = 0, resistiveLpRight = 0, lpCoeff = 0, resistiveLp2Left = 0, resistiveLp2Right = 0, resistiveLp3Left = 0, resistiveLp3Right = 0, lp2Coeff = 0;
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
    struct Frame { std::vector<float> wr, wi, qr, qi, hr, hi, omega, sigma, rad, phi; std::array<double, kZones> accel {}, lastVel {}, force {}; double knock = 0, knockDecay = 0; } frame;
    void buildFrame();
    double radScale = 1.0;   // output units of the body model's radiation (set in buildSoundboard)
    void runFrame (double& left, double& right);

    // Neighbouring strings whose modes sit back to back are stepped in one go
    struct Run { int first, count; float panLeft, panRight; bool side; int zone; };
    std::vector<Run> runs;
    std::vector<float> runGain, runPanLeft, runPanRight;   // each stepped run's gain on this sample's bridge acceleration, and its pan
    std::vector<int> runDim;                               // ...and which bridge zone (vertical, or kZones + zone sideways) it pushes
    std::vector<double> bridgeY;              // the board's same-sample motion per unit force, zone by zone (vertical then sideways)
    void buildBridgeAdmittance();
    struct Segment { int first, count; float panLeft, panRight; bool side; int zone; float topF; };
    std::vector<Segment> segmentScratch;
};

} // namespace openpiano::piano
