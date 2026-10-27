/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include <juce_gui_basics/juce_gui_basics.h>

namespace sky
{
    // Paints the dusk sky: night fading to a sunset band, faint streaks of cloud, a few stars,
    // a crescent moon and a black mountain ridge. Expensive: render once into an image and reuse it.
    void paint (juce::Graphics& g, juce::Rectangle<float> area);

    // A small crescent moon for the wordmark.
    void paintMoonIcon (juce::Graphics& g, juce::Rectangle<float> area, juce::Colour colour);
}
