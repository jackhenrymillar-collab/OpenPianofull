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

void paint (Graphics& g, Rectangle<float> area)
{
    const float w = area.getWidth(), h = area.getHeight();

    // Sky: deep blue overhead, warming to a golden horizon at the bottom
    ColourGradient skyFill (tc::skyTop, area.getTopLeft(), tc::horizon, area.getBottomLeft(), false);
    skyFill.addColour (0.55, tc::skyMid);
    skyFill.addColour (0.82, tc::skyMid.interpolatedWith (tc::horizon, 0.55f));
    g.setGradientFill (skyFill);
    g.fillRect (area);

    // The sun, centred below the bottom edge so only its upper rim shows
    const Point<float> sun (area.getCentreX(), area.getBottom() + h * 0.10f);
    const float sunR = h * 0.20f;

    // Wide atmospheric glow
    softBlob (g, sun, w * 0.75f, h * 0.85f, tc::sunGlow.withAlpha (0.55f));
    softBlob (g, sun, w * 0.42f, h * 0.55f, tc::sunGlow.withAlpha (0.45f));

    // Long, faint rays fanning upwards
    {
        Graphics::ScopedSaveState s (g);
        const int rays = 14;
        for (int i = 0; i < rays; ++i)
        {
            const float angle = MathConstants<float>::pi * (1.08f + 0.84f * (float) i / (rays - 1)); // left to right, upper half
            const float spread = 0.035f + 0.02f * (float) (i % 3);
            Path ray;
            ray.startNewSubPath (sun);
            ray.lineTo (sun + Point<float> (std::cos (angle - spread), std::sin (angle - spread)) * h * 1.3f);
            ray.lineTo (sun + Point<float> (std::cos (angle + spread), std::sin (angle + spread)) * h * 1.3f);
            ray.closeSubPath();
            ColourGradient fade (tc::sunCore.withAlpha (0.16f), sun, tc::sunCore.withAlpha (0.0f),
                                 sun.translated (0, -h * 1.0f), true);
            g.setGradientFill (fade);
            g.fillPath (ray);
        }
    }

    // Disc and bright corona
    softBlob (g, sun, sunR * 1.9f, sunR * 1.9f, tc::sunCore.withAlpha (0.65f));
    ColourGradient disc (Colours::white, sun.translated (0, -sunR * 0.3f), tc::sunCore, sun.translated (0, -sunR), true);
    g.setGradientFill (disc);
    g.fillEllipse (sun.x - sunR, sun.y - sunR, sunR * 2, sunR * 2);

    // Clouds drifting across the upper sky, lit warmly lower down
    cloud (g, { w * 0.10f, h * 0.10f }, 70.0f, 0.38f);
    cloud (g, { w * 0.88f, h * 0.07f }, 90.0f, 0.34f);
    cloud (g, { w * 0.62f, h * 0.15f }, 55.0f, 0.22f);
    cloud (g, { w * 0.03f, h * 0.55f }, 80.0f, 0.30f);
    cloud (g, { w * 0.97f, h * 0.50f }, 95.0f, 0.32f);

    // A low haze band where the light is strongest
    softBlob (g, { area.getCentreX(), area.getBottom() - h * 0.08f }, w * 0.7f, h * 0.10f, Colours::white.withAlpha (0.18f));
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
