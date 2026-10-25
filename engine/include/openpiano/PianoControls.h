/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

namespace openpiano
{

// Player-facing controls for the piano. 1.0 is the neutral setting for each.
struct PianoControls
{
    float outputGain = 1.0f;      // Linear gain on top of the engine's internal level
    float dynamicsCurve = 1.0f;   // Velocity -> hammer speed exponent (1 = linear, >1 = softer touch)
    float hammerHardness = 1.0f;  // Scales felt stiffness (brighter above 1)
    float damperStrength = 1.0f;  // Scales how quickly released notes are damped
    float actionNoise = 1.0f;     // Level of key, hammer and damper noises (0 = off)
};

} // namespace openpiano
