/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include <juce_gui_basics/juce_gui_basics.h>

// Dusk variant of the Overcast family look: a sky fading from night to a sunset band over a
// mountain ridge, smoked-glass panels, moonlit pearl knobs, serif display type, spaced small
// caps and sunset amber for anything selected or lit.
namespace theme
{
    namespace colours
    {
        // Dusk sky, sampled from Jack's reference photo: night overhead, slate and mauve
        // in the middle, a band of sunset orange and gold above a black mountain ridge
        inline const juce::Colour skyTop      { 0xff080b1c }; // night overhead
        inline const juce::Colour skyHigh     { 0xff141a3c }; // deep indigo
        inline const juce::Colour skyMid      { 0xff3b3f62 }; // slate violet
        inline const juce::Colour skyLow      { 0xff7c7488 }; // mauve grey
        inline const juce::Colour peach       { 0xffd99468 };
        inline const juce::Colour sunset      { 0xffef7424 }; // the orange band
        inline const juce::Colour gold        { 0xfff2a83a }; // just above the ridge
        inline const juce::Colour ridge       { 0xff050508 }; // mountain silhouette
        inline const juce::Colour moon        { 0xfffff6e0 };
        inline const juce::Colour sunGlow     { 0xffff9a40 }; // light from lit keys
        inline const juce::Colour amber       { 0xffff9a3c }; // selected / lit
        inline const juce::Colour ink         { 0xfff1eee8 }; // main text (light on dusk)
        inline const juce::Colour inkDim      { 0xffa7aac4 };
        inline const juce::Colour shade       { 0xff03040c }; // shadows
        inline const juce::Colour glass       { 0x660c1030 }; // smoked panel fill
        inline const juce::Colour glassBorder { 0x40ffffff };
        inline const juce::Colour pearlLight  { 0xfff4f2ef };
        inline const juce::Colour pearlDark   { 0xff8e93ad };
    }

    // Cormorant Garamond: wordmark
    juce::Font serif (float height, float tracking = 0.0f);

    // Jost: labels, values and buttons
    juce::Font sans (float height, float tracking = 0.0f);

    // Spaced small caps label, e.g. "VOLUME"
    juce::Font label (float height = 12.0f);

    // A soft glow around a rounded rectangle (for selected controls).
    void drawGlow (juce::Graphics& g, juce::Rectangle<float> r, float cornerRadius, juce::Colour colour, float size = 8.0f);

    // A smoked-glass panel: translucent dark with a faint hairline edge.
    void drawGlassPanel (juce::Graphics& g, juce::Rectangle<float> r, float cornerRadius, float opacity = 1.0f);
}
