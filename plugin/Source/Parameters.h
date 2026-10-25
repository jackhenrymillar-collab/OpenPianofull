/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include <juce_audio_processors/juce_audio_processors.h>

namespace params
{
    inline constexpr auto volume = "volume";
    inline constexpr auto dynamics = "dynamics";
    inline constexpr auto hardness = "hardness";
    inline constexpr auto damper = "damper";
    inline constexpr auto actionNoise = "actionNoise";
    inline constexpr auto age = "age";
    inline constexpr auto sustainLatch = "sustainLatch";
    inline constexpr auto velocitySens = "velocitySens";

    juce::AudioProcessorValueTreeState::ParameterLayout createLayout();

    struct FactoryPreset
    {
        const char* name;
        float volumeDb, dynamics, hardness, damper;
    };

    const std::vector<FactoryPreset>& factoryPresets();
}
