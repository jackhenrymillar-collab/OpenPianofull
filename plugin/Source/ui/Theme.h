/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include <juce_gui_basics/juce_gui_basics.h>

// Daylight variant of the Overcast family look: a bright sky with the sun rising from below,
// frosted-glass panels, pearl knobs, serif display type, spaced small caps and warm
// sunlight amber for anything selected or lit.
namespace theme
{
    namespace colours
    {
        inline const juce::Colour skyTop      { 0xff1846c4 }; // deep summer blue overhead
        inline const juce::Colour skyMid      { 0xff2f74e0 };
        inline const juce::Colour horizon     { 0xff8cc4f6 }; // paler near the field
        inline const juce::Colour field       { 0xfff2d01e }; // rapeseed yellow
        inline const juce::Colour fieldDeep   { 0xffb89a0c };
        inline const juce::Colour stem        { 0xff5f8f2a };
        inline const juce::Colour sunCore     { 0xfffff8e2 };
        inline const juce::Colour sunGlow     { 0xffffc46a };
        inline const juce::Colour amber       { 0xffef9b2c }; // selected / lit
        inline const juce::Colour ink         { 0xff1c2a52 }; // main text
        inline const juce::Colour inkDim      { 0xff4f5f8c };
        inline const juce::Colour glass       { 0x66ffffff }; // frosted panel fill
        inline const juce::Colour glassBorder { 0xccffffff };
        inline const juce::Colour pearlLight  { 0xffffffff };
        inline const juce::Colour pearlDark   { 0xffcfd6e6 };
    }

    // Cormorant Garamond: wordmark
    juce::Font serif (float height, float tracking = 0.0f);

    // Jost: labels, values and buttons
    juce::Font sans (float height, float tracking = 0.0f);

    // Spaced small caps label, e.g. "VOLUME"
    juce::Font label (float height = 12.0f);

    // A soft glow around a rounded rectangle (for selected controls).
    void drawGlow (juce::Graphics& g, juce::Rectangle<float> r, float cornerRadius, juce::Colour colour, float size = 8.0f);

    // A frosted-glass panel: translucent white with a bright hairline edge.
    void drawGlassPanel (juce::Graphics& g, juce::Rectangle<float> r, float cornerRadius, float opacity = 1.0f);
}
