/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#include "Theme.h"
#include "BinaryData.h"

namespace theme
{

static juce::Typeface::Ptr loadTypeface (const void* data, size_t size)
{
    return juce::Typeface::createSystemTypefaceFor (data, size);
}

static juce::Typeface::Ptr serifTypeface()
{
    static auto tf = loadTypeface (OpenPianoAssets::CormorantGaramond_ttf, OpenPianoAssets::CormorantGaramond_ttfSize);
    return tf;
}

static juce::Typeface::Ptr sansTypeface()
{
    static auto tf = loadTypeface (OpenPianoAssets::Jost_ttf, OpenPianoAssets::Jost_ttfSize);
    return tf;
}

juce::Font serif (float height, float tracking)
{
    return juce::Font (juce::FontOptions (serifTypeface()).withHeight (height).withKerningFactor (tracking));
}

juce::Font sans (float height, float tracking)
{
    return juce::Font (juce::FontOptions (sansTypeface()).withHeight (height).withKerningFactor (tracking));
}

juce::Font label (float height)
{
    return sans (height, 0.22f);
}

void drawGlow (juce::Graphics& g, juce::Rectangle<float> r, float cornerRadius, juce::Colour colour, float size)
{
    const int steps = 8;
    for (int i = steps; i >= 1; --i)
    {
        const float grow = size * (float) i / steps;
        const float alpha = 0.10f * (1.0f - (float) (i - 1) / steps);
        g.setColour (colour.withAlpha (alpha));
        g.drawRoundedRectangle (r.expanded (grow * 0.5f), cornerRadius + grow * 0.5f, grow);
    }
}

void drawGlassPanel (juce::Graphics& g, juce::Rectangle<float> r, float cornerRadius, float opacity)
{
    // Soft shadow underneath
    g.setColour (colours::shade.withAlpha (0.25f * opacity));
    g.fillRoundedRectangle (r.translated (0, 3.0f).expanded (1.0f), cornerRadius + 1.0f);

    // Smoked glass: a little lighter at the top, where the sky is darkest
    juce::ColourGradient fill (juce::Colour (0xff1c2148).withAlpha (0.42f * opacity), r.getTopLeft(),
                               juce::Colour (0xff0a0c1e).withAlpha (0.38f * opacity), r.getBottomLeft(), false);
    g.setGradientFill (fill);
    g.fillRoundedRectangle (r, cornerRadius);

    // Hairline edge, catching a little warm light along the bottom
    juce::ColourGradient edge (juce::Colours::white.withAlpha (0.28f * opacity), r.getTopLeft(),
                               colours::sunset.withAlpha (0.35f * opacity), r.getBottomLeft(), false);
    g.setGradientFill (edge);
    g.drawRoundedRectangle (r.reduced (0.5f), cornerRadius, 1.0f);
}

} // namespace theme
