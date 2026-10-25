/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include <juce_gui_basics/juce_gui_basics.h>

// Visual language shared with the Overcast plugin family: a night sky, serif display type,
// spaced small caps, glassy navy controls and warm gold for anything selected.
namespace theme
{
    namespace colours
    {
        inline const juce::Colour skyTop      { 0xff0b1030 };
        inline const juce::Colour skyMid      { 0xff1a2462 };
        inline const juce::Colour skyLow      { 0xff2c3780 };
        inline const juce::Colour moon        { 0xfff4ecc8 };
        inline const juce::Colour moonShade   { 0xffd9cc98 };
        inline const juce::Colour gold        { 0xffe9c46a };
        inline const juce::Colour text        { 0xfff3f1ff };
        inline const juce::Colour textDim     { 0xffa9add6 };
        inline const juce::Colour panel       { 0x8c141a46 }; // translucent navy
        inline const juce::Colour panelBorder { 0x33ffffff };
        inline const juce::Colour knobLight   { 0xff3d4ba4 };
        inline const juce::Colour knobDark    { 0xff161d58 };
    }

    // Cormorant Garamond: wordmark and titles
    juce::Font serif (float height, float tracking = 0.0f);

    // Jost: labels, values and buttons
    juce::Font sans (float height, float tracking = 0.0f);

    // Spaced small caps label, e.g. "DECAY"
    juce::Font label (float height = 12.0f);

    // Draws a soft glow around a rounded rectangle (for selected controls).
    void drawGlow (juce::Graphics& g, juce::Rectangle<float> r, float cornerRadius, juce::Colour colour, float size = 8.0f);
}
