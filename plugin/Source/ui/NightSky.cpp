/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#include "NightSky.h"
#include "Theme.h"

namespace nightsky
{

using namespace juce;

static void softBlob (Graphics& g, Point<float> centre, float rx, float ry, Colour colour)
{
    // A radial gradient stretched into an ellipse: a cheap, soft cloud puff
    Graphics::ScopedSaveState s (g);
    g.addTransform (AffineTransform::scale (1.0f, ry / rx, centre.x, centre.y));
    ColourGradient grad (colour, centre, colour.withAlpha (0.0f), centre.translated (rx, 0), true);
    g.setGradientFill (grad);
    g.fillEllipse (centre.x - rx, centre.y - rx, rx * 2, rx * 2);
}

static void paintStars (Graphics& g, Rectangle<float> area)
{
    Random rng (20261025); // fixed seed: the same sky every time
    const int count = (int) (area.getWidth() * area.getHeight() / 2600.0f);
    for (int i = 0; i < count; ++i)
    {
        const float x = area.getX() + rng.nextFloat() * area.getWidth();
        const float y = area.getY() + rng.nextFloat() * area.getHeight() * 0.85f;
        const float size = 0.6f + rng.nextFloat() * rng.nextFloat() * 1.8f;
        const float alpha = 0.15f + rng.nextFloat() * 0.55f;
        g.setColour (Colours::white.withAlpha (alpha));
        g.fillEllipse (x - size * 0.5f, y - size * 0.5f, size, size);

        if (size > 1.8f) // a faint halo on the brightest ones
            softBlob (g, { x, y }, size * 3.0f, size * 3.0f, Colours::white.withAlpha (alpha * 0.25f));
    }
}

static void paintMoon (Graphics& g, Point<float> c, float r)
{
    // Halo
    softBlob (g, c, r * 2.0f, r * 2.0f, theme::colours::moon.withAlpha (0.16f));
    softBlob (g, c, r * 1.35f, r * 1.35f, theme::colours::moon.withAlpha (0.22f));

    // Disc, lit from the upper left
    ColourGradient disc (theme::colours::moon.brighter (0.05f), c.translated (-r * 0.4f, -r * 0.45f),
                         theme::colours::moonShade, c.translated (r * 0.7f, r * 0.8f), true);
    g.setGradientFill (disc);
    g.fillEllipse (c.x - r, c.y - r, r * 2, r * 2);

    // Craters
    struct Crater { float x, y, size; };
    const Crater craters[] = { { -0.32f, -0.22f, 0.17f }, { 0.36f, 0.16f, 0.22f }, { 0.36f, -0.48f, 0.08f }, { -0.06f, 0.48f, 0.10f } };
    for (auto& cr : craters)
    {
        const auto cc = c.translated (cr.x * r, cr.y * r);
        const float cs = cr.size * r;
        g.setColour (theme::colours::moonShade.darker (0.12f).withAlpha (0.75f));
        g.fillEllipse (cc.x - cs, cc.y - cs, cs * 2, cs * 2);
    }
}

void paint (Graphics& g, Rectangle<float> area, const Layout& layout)
{
    // Sky
    ColourGradient sky (theme::colours::skyTop, area.getTopLeft(), theme::colours::skyLow, area.getBottomLeft(), false);
    sky.addColour (0.55, theme::colours::skyMid);
    g.setGradientFill (sky);
    g.fillRect (area);

    // A faint lighter wash behind the moon
    softBlob (g, layout.moonCentre, area.getWidth() * 0.45f, area.getHeight() * 0.5f, Colour (0x223a4cb0));

    paintStars (g, area);
    paintMoon (g, layout.moonCentre, layout.moonRadius);

    // A wisp of cloud drifting past the moon
    const auto mc = layout.moonCentre;
    const float mr = layout.moonRadius;
    softBlob (g, mc.translated (-mr * 1.2f, mr * 0.6f), mr * 1.3f, mr * 0.35f, Colour (0x30c8cdf5));
    softBlob (g, mc.translated (-mr * 0.4f, mr * 0.75f), mr * 0.9f, mr * 0.25f, Colour (0x28d6daf8));

    // Cloud bank along the bottom
    Random rng (7);
    const float base = area.getBottom() - area.getHeight() * 0.20f;
    for (int i = 0; i < 26; ++i)
    {
        const float x = area.getX() + area.getWidth() * (i / 25.0f) + (rng.nextFloat() - 0.5f) * 60.0f;
        const float y = base + rng.nextFloat() * area.getHeight() * 0.18f;
        const float rx = 90.0f + rng.nextFloat() * 120.0f;
        const float ry = 28.0f + rng.nextFloat() * 32.0f;
        softBlob (g, { x, y }, rx, ry, Colour (0xffb9bff0).withAlpha (0.10f + rng.nextFloat() * 0.10f));
    }

    // Darken the very bottom slightly so the keyboard sits on solid ground
    ColourGradient floor (Colours::transparentBlack, 0, area.getBottom() - 140.0f,
                          Colour (0x660a0e2a), 0, area.getBottom(), false);
    g.setGradientFill (floor);
    g.fillRect (area.withTop (area.getBottom() - 140.0f));
}

void paintCrescent (Graphics& g, Rectangle<float> area, Colour colour)
{
    const auto d = jmin (area.getWidth(), area.getHeight());
    const auto disc = area.withSizeKeepingCentre (d, d);
    // Disc and bite in one even-odd path: clipping to it and filling the disc leaves the
    // disc minus the bite.
    Path p;
    p.addEllipse (disc);
    p.addEllipse (disc.translated (d * 0.32f, -d * 0.12f));
    p.setUsingNonZeroWinding (false);
    Graphics::ScopedSaveState s (g);
    g.reduceClipRegion (p);
    g.setColour (colour);
    g.fillEllipse (disc);
}

} // namespace nightsky
