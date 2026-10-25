/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include <juce_gui_basics/juce_gui_basics.h>

namespace nightsky
{
    struct Layout
    {
        juce::Point<float> moonCentre;
        float moonRadius;
    };

    // Paints the sky gradient, stars, moon and clouds. Expensive: render once into an
    // image and reuse it.
    void paint (juce::Graphics& g, juce::Rectangle<float> area, const Layout& layout);

    // Draws a small crescent moon icon (for the wordmark).
    void paintCrescent (juce::Graphics& g, juce::Rectangle<float> area, juce::Colour colour);
}
