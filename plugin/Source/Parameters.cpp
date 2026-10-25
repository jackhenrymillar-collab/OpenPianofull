/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#include "Parameters.h"

namespace params
{

using namespace juce;

static NormalisableRange<float> scaleRange (float lo, float hi)
{
    // Logarithmic, centred on 1.0 (100%)
    NormalisableRange<float> r (lo, hi);
    r.setSkewForCentre (1.0f);
    return r;
}

AudioProcessorValueTreeState::ParameterLayout createLayout()
{
    AudioProcessorValueTreeState::ParameterLayout layout;

    layout.add (std::make_unique<AudioParameterFloat> (
        ParameterID { volume, 1 }, "Volume", NormalisableRange<float> (-36.0f, 12.0f, 0.1f), 0.0f,
        AudioParameterFloatAttributes()
            .withLabel ("dB")
            .withStringFromValueFunction ([] (float v, int) { return String (v, 1) + " dB"; })));

    layout.add (std::make_unique<AudioParameterFloat> (
        ParameterID { dynamics, 1 }, "Dynamics", scaleRange (0.5f, 2.0f), 1.0f,
        AudioParameterFloatAttributes()
            .withStringFromValueFunction ([] (float v, int)
            {
                // Above 1, soft playing gets softer. Below 1, the range is compressed.
                if (std::abs (v - 1.0f) < 0.005f)
                    return String ("Linear");
                return String (v > 1.0f ? "Wide " : "Narrow ") + String (v, 2);
            })));

    auto percent = [] (float v, int) { return String (roundToInt (v * 100.0f)) + "%"; };

    layout.add (std::make_unique<AudioParameterFloat> (
        ParameterID { hardness, 1 }, "Hammer Hardness", scaleRange (0.25f, 4.0f), 1.0f,
        AudioParameterFloatAttributes().withStringFromValueFunction (percent)));

    layout.add (std::make_unique<AudioParameterFloat> (
        ParameterID { damper, 1 }, "Damper", scaleRange (0.25f, 4.0f), 1.0f,
        AudioParameterFloatAttributes().withStringFromValueFunction (percent)));

    layout.add (std::make_unique<AudioParameterBool> (ParameterID { sustainLatch, 1 }, "Sustain Latch", false));

    return layout;
}

const std::vector<FactoryPreset>& factoryPresets()
{
    static const std::vector<FactoryPreset> presets {
        { "Open Grand",     0.0f, 1.0f, 1.0f,  1.0f },
        { "Soft Felt",      2.0f, 1.3f, 0.45f, 1.0f },
        { "Bright Hammers", -3.0f, 1.0f, 2.2f, 1.0f },
        { "Gentle Touch",   0.0f, 1.6f, 0.8f,  1.0f },
        { "Quick Dampers",  0.0f, 1.0f, 1.0f,  3.0f },
    };
    return presets;
}

} // namespace params
