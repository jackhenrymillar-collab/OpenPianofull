/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#include "Sky.h"
#include "Theme.h"

namespace sky
{

using namespace juce;
namespace tc = theme::colours;

static void softBlob (Graphics& g, Point<float> centre, float rx, float ry, Colour colour)
{
    // A radial gradient stretched into an ellipse
    Graphics::ScopedSaveState s (g);
    g.addTransform (AffineTransform::scale (1.0f, ry / rx, centre.x, centre.y));
    ColourGradient grad (colour, centre, colour.withAlpha (0.0f), centre.translated (rx, 0), true);
    g.setGradientFill (grad);
    g.fillEllipse (centre.x - rx, centre.y - rx, rx * 2, rx * 2);
}

static void crescent (Graphics& g, Point<float> c, float r, Colour colour)
{
    // The lit limb on the lower right, the shadowed part offset up and to the left
    Path disc;
    disc.addEllipse (c.x - r, c.y - r, r * 2, r * 2);
    Path p (disc);
    p.addEllipse (c.x - r * 1.32f, c.y - r * 1.38f, r * 2.05f, r * 2.05f);
    p.setUsingNonZeroWinding (false);
    Graphics::ScopedSaveState s (g);
    g.reduceClipRegion (disc);
    g.setColour (colour);
    g.fillPath (p);
}

static void ridge (Graphics& g, Rectangle<float> area, float baseY)
{
    // A low, jagged mountain line, black against the sunset, with one small peak off centre
    const float w = area.getWidth();
    Random rng (4021);
    Path p;
    p.startNewSubPath (area.getX(), area.getBottom());
    float y = baseY + 2.0f;
    for (float x = 0.0f; x <= w; x += 6.0f)
    {
        y += (rng.nextFloat() - 0.5f) * 2.4f;
        y = jlimit (baseY - 7.0f, baseY + 5.0f, y + (baseY - y) * 0.08f);
        float bump = 0.0f;
        const auto peak = [&] (float px, float width, float height)
        {
            const float d = (x - px) / width;
            if (std::abs (d) < 1.0f) bump -= height * (1.0f - std::abs (d)) * (0.75f + 0.5f * rng.nextFloat());
        };
        peak (w * 0.50f, 26.0f, 10.0f);
        peak (w * 0.365f, 14.0f, 6.0f);
        peak (w * 0.24f, 40.0f, 5.0f);
        peak (w * 0.82f, 60.0f, 6.0f);
        p.lineTo (area.getX() + x, y + bump);
    }
    p.lineTo (area.getRight(), area.getBottom());
    p.closeSubPath();
    g.setColour (tc::ridge);
    g.fillPath (p);
}

void paint (Graphics& g, Rectangle<float> area)
{
    const float w = area.getWidth(), h = area.getHeight();
    const float horizonY = area.getY() + h * 0.905f;

    // Sky: night overhead, through indigo, slate and mauve to a sunset band at the ridge
    ColourGradient skyFill (tc::skyTop, area.getTopLeft(), tc::gold, { area.getX(), horizonY }, false);
    skyFill.addColour (0.22, tc::skyHigh);
    skyFill.addColour (0.50, tc::skyMid);
    skyFill.addColour (0.72, tc::skyLow);
    skyFill.addColour (0.84, tc::peach);
    skyFill.addColour (0.93, tc::sunset);
    g.setGradientFill (skyFill);
    g.fillRect (area);

    // The sunset is a little brighter and more golden to the right, as in the photo
    softBlob (g, { w * 0.86f, horizonY - 6.0f }, w * 0.45f, 46.0f, tc::gold.withAlpha (0.55f));
    softBlob (g, { w * 0.25f, horizonY - 4.0f }, w * 0.40f, 22.0f, Colour (0xffe0501c).withAlpha (0.45f));

    // Faint streaks of high cloud, stretched sideways and slightly tilted
    {
        Random rng (911);
        for (int i = 0; i < 26; ++i)
        {
            const float y = h * (0.12f + 0.62f * rng.nextFloat());
            const float x = w * rng.nextFloat();
            Graphics::ScopedSaveState s (g);
            g.addTransform (AffineTransform::rotation (-0.06f, x, y));
            const bool low = y > h * 0.5f;
            const auto tint = low ? Colour (0xffb08a8a) : Colour (0xff2c3260);
            softBlob (g, { x, y }, w * (0.18f + 0.25f * rng.nextFloat()), 4.0f + 6.0f * rng.nextFloat(),
                      tint.withAlpha (0.10f + 0.10f * rng.nextFloat()));
        }
    }

    // A thin dark haze line just above the ridge
    softBlob (g, { w * 0.30f, horizonY - 3.0f }, w * 0.35f, 2.5f, Colour (0xff5a2a20).withAlpha (0.5f));

    // A few faint stars in the darkest part
    {
        Random rng (5);
        for (int i = 0; i < 60; ++i)
        {
            const float y = h * 0.30f * rng.nextFloat() * rng.nextFloat() * 1.6f;
            const float r = 0.4f + 0.7f * rng.nextFloat();
            g.setColour (Colours::white.withAlpha ((0.15f + 0.45f * rng.nextFloat()) * (1.0f - y / (h * 0.5f))));
            g.fillEllipse (w * rng.nextFloat() - r, y - r, r * 2, r * 2);
        }
    }

    // The crescent moon, high up between the wordmark and the mode buttons
    const Point<float> moonCentre (w * 0.585f, h * 0.058f);
    softBlob (g, moonCentre, 46.0f, 46.0f, tc::moon.withAlpha (0.10f));
    softBlob (g, moonCentre, 18.0f, 18.0f, tc::moon.withAlpha (0.16f));
    crescent (g, moonCentre, 9.0f, tc::moon);

    ridge (g, area, horizonY);
}

void paintMoonIcon (Graphics& g, Rectangle<float> area, Colour colour)
{
    const float r = jmin (area.getWidth(), area.getHeight()) * 0.42f;
    crescent (g, area.getCentre(), r, colour);
}

} // namespace sky
