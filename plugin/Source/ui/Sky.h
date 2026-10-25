/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include <juce_gui_basics/juce_gui_basics.h>

namespace sky
{
    // Paints the daylight sky: gradient, drifting clouds, and the sun peeking up from
    // below the bottom edge. Expensive: render once into an image and reuse it.
    void paint (juce::Graphics& g, juce::Rectangle<float> area);

    // A small sun icon for the wordmark.
    void paintSunIcon (juce::Graphics& g, juce::Rectangle<float> area, juce::Colour colour);
}
