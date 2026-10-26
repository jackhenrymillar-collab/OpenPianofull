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

static void cloud (Graphics& g, Point<float> c, float size, float alpha)
{
    // A cumulus made of overlapping puffs, flatter at the base
    Random rng ((int64) (c.x * 31 + c.y * 17));
    for (int i = 0; i < 9; ++i)
    {
        const float dx = (rng.nextFloat() - 0.5f) * size * 1.6f;
        const float dy = -rng.nextFloat() * size * 0.35f;
        const float r = size * (0.35f + rng.nextFloat() * 0.3f);
        softBlob (g, c.translated (dx, dy), r, r * 0.75f, Colours::white.withAlpha (alpha));
    }
}

static void cumulus (Graphics& g, Point<float> c, float size, float alpha)
{
    // A bright, crisp-edged summer cumulus: many small dense puffs over a few big ones,
    // with a faint blue-grey shade along the flat base
    Random rng ((int64) (c.x * 31 + c.y * 17));
    for (int i = 0; i < 7; ++i)
    {
        const float dx = (rng.nextFloat() - 0.5f) * size * 1.7f;
        const float r = size * (0.32f + rng.nextFloat() * 0.25f);
        softBlob (g, c.translated (dx, size * 0.18f), r, r * 0.55f, Colour (0xffc9d6ea).withAlpha (alpha * 0.55f));
    }
    for (int i = 0; i < 16; ++i)
    {
        const float dx = (rng.nextFloat() - 0.5f) * size * 1.6f;
        const float dy = -rng.nextFloat() * size * 0.55f;
        const float r = size * (0.22f + rng.nextFloat() * 0.25f);
        softBlob (g, c.translated (dx, dy), r, r * 0.9f, Colours::white.withAlpha (alpha));
        softBlob (g, c.translated (dx, dy), r * 0.6f, r * 0.55f, Colours::white.withAlpha (alpha * 0.8f));
    }
}

static void wisp (Graphics& g, Point<float> c, float size, float alpha)
{
    Random rng ((int64) (c.x * 13 + c.y * 7));
    for (int i = 0; i < 6; ++i)
        softBlob (g, c.translated ((rng.nextFloat() - 0.5f) * size * 1.8f, (rng.nextFloat() - 0.5f) * size * 0.5f),
                  size * (0.3f + 0.2f * rng.nextFloat()), size * 0.18f, Colours::white.withAlpha (alpha));
}

static void field (Graphics& g, Rectangle<float> area)
{
    // A field of rapeseed in flower, seen from inside it: a soft yellow mass at the bottom,
    // then stems and clusters of four-petalled flowers rising against the sky
    const float w = area.getWidth(), h = area.getHeight();
    const float top = area.getBottom() - h * 0.17f;
    Random rng (20240607);

    // Distant, out-of-focus yellow
    softBlob (g, { w * 0.5f, area.getBottom() - h * 0.05f }, w * 0.85f, h * 0.12f, tc::field.withAlpha (0.85f));
    {
        ColourGradient base (tc::field.withAlpha (0.0f), 0.0f, top + h * 0.06f, tc::field, 0.0f, area.getBottom() - h * 0.08f, false);
        base.addColour (0.7, tc::field.interpolatedWith (tc::fieldDeep, 0.25f));
        g.setGradientFill (base);
        g.fillRect (area.withTop (top + h * 0.06f));
    }

    // Stems swaying a little, each topped with a loose cluster of flowers
    for (int i = 0; i < 90; ++i)
    {
        const float x = rng.nextFloat() * w;
        const float tipY = top + rng.nextFloat() * h * 0.10f;
        const float lean = (rng.nextFloat() - 0.5f) * 30.0f;
        Path stem;
        stem.startNewSubPath (x, area.getBottom());
        stem.quadraticTo (x + lean * 0.3f, (tipY + area.getBottom()) * 0.5f, x + lean, tipY);
        g.setColour (tc::stem.withAlpha (0.55f + 0.3f * rng.nextFloat()));
        g.strokePath (stem, PathStrokeType (1.0f + rng.nextFloat()));

        const int flowers = 4 + rng.nextInt (6);
        for (int f = 0; f < flowers; ++f)
        {
            const auto p = Point<float> (x + lean + (rng.nextFloat() - 0.5f) * 22.0f, tipY + rng.nextFloat() * 30.0f - 6.0f);
            const float r = 2.6f + rng.nextFloat() * 3.4f;
            const auto petal = tc::field.interpolatedWith (Colour (0xfffff27a), rng.nextFloat() * 0.6f);
            for (int k = 0; k < 4; ++k)
            {
                const float a = MathConstants<float>::halfPi * (float) k + rng.nextFloat() * 0.4f;
                g.setColour (petal);
                g.fillEllipse (Rectangle<float> (r * 1.1f, r * 1.1f).withCentre (p.getPointOnCircumference (r * 0.5f, a)));
            }
            g.setColour (tc::fieldDeep.withAlpha (0.7f));
            g.fillEllipse (Rectangle<float> (r * 0.5f, r * 0.5f).withCentre (p));
        }
    }
}

void paint (Graphics& g, Rectangle<float> area)
{
    const float w = area.getWidth(), h = area.getHeight();

    // Sky: deep saturated blue overhead, paler towards the horizon
    ColourGradient skyFill (tc::skyTop, area.getTopLeft(), tc::horizon, area.getBottomLeft(), false);
    skyFill.addColour (0.45, tc::skyMid);
    g.setGradientFill (skyFill);
    g.fillRect (area);

    // The sun, high in the sky just above the keyboard, with a bright starburst
    const Point<float> sun (w * 0.47f, h * 0.075f);
    const float sunR = 15.0f;
    softBlob (g, sun, w * 0.40f, h * 0.40f, Colour (0xff9ccaff).withAlpha (0.45f)); // blue-white haze
    softBlob (g, sun, 120.0f, 120.0f, Colours::white.withAlpha (0.45f));

    // Sparkle: thin tapering rays of different lengths
    {
        Random rng (77);
        const int rays = 28;
        for (int i = 0; i < rays; ++i)
        {
            const float angle = MathConstants<float>::twoPi * ((float) i + 0.3f * rng.nextFloat()) / rays;
            const float len = (i % 2 == 0 ? 150.0f : 85.0f) * (0.7f + 0.6f * rng.nextFloat());
            const float half = 0.022f + 0.01f * rng.nextFloat();
            Path ray;
            ray.startNewSubPath (sun.getPointOnCircumference (sunR * 0.6f, angle - half * 6.0f));
            ray.lineTo (sun.getPointOnCircumference (len, angle));
            ray.lineTo (sun.getPointOnCircumference (sunR * 0.6f, angle + half * 6.0f));
            ray.closeSubPath();
            ColourGradient fade (Colours::white.withAlpha (0.85f), sun, Colours::white.withAlpha (0.0f),
                                 sun.getPointOnCircumference (len, angle), true);
            g.setGradientFill (fade);
            g.fillPath (ray);
        }
    }
    softBlob (g, sun, sunR * 3.2f, sunR * 3.2f, Colours::white.withAlpha (0.8f));
    g.setColour (Colours::white);
    g.fillEllipse (sun.x - sunR, sun.y - sunR, sunR * 2, sunR * 2);

    // Clouds: a big cumulus on the left, a smaller one on the right, a few wisps up high
    cumulus (g, { w * 0.09f, h * 0.36f }, 120.0f, 0.55f);
    cumulus (g, { w * 0.97f, h * 0.48f }, 70.0f, 0.5f);
    wisp (g, { w * 0.86f, h * 0.05f }, 60.0f, 0.35f);
    wisp (g, { w * 0.70f, h * 0.10f }, 30.0f, 0.3f);

    field (g, area);
}

void paintSunIcon (Graphics& g, Rectangle<float> area, Colour colour)
{
    const auto c = area.getCentre();
    const float r = jmin (area.getWidth(), area.getHeight()) * 0.5f;
    g.setColour (colour);
    g.fillEllipse (Rectangle<float> (r * 0.9f, r * 0.9f).withCentre (c));
    for (int i = 0; i < 8; ++i)
    {
        const float a = MathConstants<float>::twoPi * (float) i / 8.0f;
        const auto p1 = c.getPointOnCircumference (r * 0.68f, a);
        const auto p2 = c.getPointOnCircumference (r * 1.0f, a);
        g.drawLine ({ p1, p2 }, 1.6f);
    }
}

} // namespace sky
