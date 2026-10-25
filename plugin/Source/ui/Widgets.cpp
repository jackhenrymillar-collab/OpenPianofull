/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#include "Widgets.h"
#include "Theme.h"

namespace widgets
{

using namespace juce;
namespace tc = theme::colours;

//==============================================================================
NightLookAndFeel::NightLookAndFeel()
{
    setColour (Slider::textBoxTextColourId, tc::text);
    setColour (TextButton::textColourOffId, tc::text);
    setColour (TextButton::textColourOnId, tc::text);
    setColour (TooltipWindow::backgroundColourId, tc::skyTop);
    setColour (TooltipWindow::textColourId, tc::text);
    setColour (TooltipWindow::outlineColourId, tc::panelBorder);
}

void NightLookAndFeel::drawRotarySlider (Graphics& g, int x, int y, int width, int height, float pos,
                                         float startAngle, float endAngle, Slider& slider)
{
    const auto bounds = Rectangle<float> ((float) x, (float) y, (float) width, (float) height).reduced (6.0f);
    const float r = jmin (bounds.getWidth(), bounds.getHeight()) * 0.5f;
    const auto c = bounds.getCentre();
    const auto disc = Rectangle<float> (r * 2, r * 2).withCentre (c);

    // Drop shadow
    for (int i = 6; i >= 1; --i)
    {
        g.setColour (Colours::black.withAlpha (0.05f));
        g.fillEllipse (disc.expanded ((float) i).translated (0, (float) i * 0.8f));
    }

    // Sphere: lit from above
    ColourGradient sphere (tc::knobLight, c.translated (-r * 0.25f, -r * 0.55f),
                           tc::knobDark, c.translated (r * 0.45f, r * 0.9f), true);
    g.setGradientFill (sphere);
    g.fillEllipse (disc);

    // Specular sheen on the upper half
    ColourGradient sheen (Colours::white.withAlpha (0.10f), c.translated (0, -r),
                          Colours::white.withAlpha (0.0f), c, false);
    g.setGradientFill (sheen);
    g.fillEllipse (disc.reduced (r * 0.08f).withTrimmedBottom (r * 0.7f));

    // Rim
    g.setColour (tc::gold.withAlpha (slider.isMouseOverOrDragging() ? 0.45f : 0.22f));
    g.drawEllipse (disc, 1.0f);

    // Glowing pointer dot
    const float angle = startAngle + pos * (endAngle - startAngle);
    const auto dot = c.getPointOnCircumference (r * 0.78f, angle);
    const float dr = jmax (3.0f, r * 0.075f);
    ColourGradient glow (tc::gold.withAlpha (0.55f), dot, tc::gold.withAlpha (0.0f), dot.translated (dr * 4.0f, 0), true);
    g.setGradientFill (glow);
    g.fillEllipse (Rectangle<float> (dr * 8, dr * 8).withCentre (dot));
    g.setColour (Colour (0xfffff6dc));
    g.fillEllipse (Rectangle<float> (dr * 2, dr * 2).withCentre (dot));
}

void NightLookAndFeel::drawButtonBackground (Graphics& g, Button& button, const Colour&, bool isHighlighted, bool isDown)
{
    // Buttons keep a margin inside their bounds so the selected glow isn't clipped
    auto r = button.getLocalBounds().toFloat().reduced (glowMargin);
    const bool pill = button.getProperties()[shapeKey].toString() == "pill";
    const float corner = pill ? r.getHeight() * 0.5f : 8.0f;
    const bool on = button.getToggleState();
    const bool enabled = button.isEnabled();

    if (on)
        theme::drawGlow (g, r, corner, tc::gold, glowMargin * 2.0f - 2.0f);

    auto fill = on ? Colour (0xb02a2f5e) : tc::panel;
    if (enabled && isHighlighted) fill = fill.brighter (0.15f);
    if (enabled && isDown) fill = fill.brighter (0.25f);
    g.setColour (fill.withMultipliedAlpha (enabled ? 1.0f : 0.6f));
    g.fillRoundedRectangle (r, corner);

    g.setColour (on ? tc::gold.withAlpha (0.9f) : tc::panelBorder.withMultipliedAlpha (enabled ? 1.0f : 0.6f));
    g.drawRoundedRectangle (r, corner, on ? 1.4f : 1.0f);
}

Font NightLookAndFeel::getTextButtonFont (TextButton& b, int height)
{
    const bool pill = b.getProperties()[shapeKey].toString() == "pill";
    return pill ? theme::label (jmin (13.0f, height * 0.36f)) : theme::sans (jmin (15.0f, height * 0.36f));
}

void NightLookAndFeel::drawButtonText (Graphics& g, TextButton& button, bool, bool)
{
    const auto badge = button.getProperties()[badgeKey].toString();
    auto area = button.getLocalBounds().toFloat();
    const float alpha = button.isEnabled() ? 1.0f : 0.45f;
    const bool pill = button.getProperties()[shapeKey].toString() == "pill";
    const auto text = pill ? button.getButtonText().toUpperCase() : button.getButtonText();

    g.setColour ((button.getToggleState() && pill ? tc::gold : tc::text).withMultipliedAlpha (alpha));
    g.setFont (getTextButtonFont (button, button.getHeight()));

    // A/B button: the active letter is gold
    if (button.getProperties().contains (activeSlotKey))
    {
        const bool b = (int) button.getProperties()[activeSlotKey] == 1;
        AttributedString as;
        as.setJustification (Justification::centred);
        const auto f = theme::sans (15.0f, 0.1f);
        as.append ("A", f, b ? tc::textDim : tc::gold);
        as.append ("  /  ", f, tc::textDim);
        as.append ("B", f, b ? tc::gold : tc::textDim);
        as.draw (g, area);
        return;
    }

    if (badge.isEmpty())
    {
        g.drawText (text, area, Justification::centred, false);
        return;
    }

    auto top = area.removeFromTop (area.getHeight() * 0.62f);
    g.drawText (text, top.withTrimmedTop (4.0f), Justification::centredBottom, false);
    g.setColour (tc::gold.withAlpha (0.75f));
    g.setFont (theme::label (9.5f));
    g.drawText (badge.toUpperCase(), area.withTrimmedBottom (4.0f), Justification::centredTop, false);
}

//==============================================================================
Knob::Knob (AudioProcessorValueTreeState& state, const String& paramId, const String& text)
    : labelText (text.toUpperCase()),
      parameter (state.getParameter (paramId)),
      attachment (state, paramId, slider)
{
    slider.setSliderStyle (Slider::RotaryVerticalDrag);
    slider.setTextBoxStyle (Slider::NoTextBox, false, 0, 0);
    slider.setRotaryParameters (MathConstants<float>::pi * 1.25f, MathConstants<float>::pi * 2.75f, true);
    slider.setDoubleClickReturnValue (true, (double) parameter->convertFrom0to1 (parameter->getDefaultValue()));
    slider.setMouseDragSensitivity (220);
    slider.onValueChange = [this] { repaint(); };
    slider.setTooltip (parameter->getName (64) + " (double-click to reset)");
    addAndMakeVisible (slider);
}

void Knob::resized()
{
    auto r = getLocalBounds();
    const int size = jmin (r.getWidth(), r.getHeight() - 52);
    slider.setBounds (r.removeFromTop (size).withSizeKeepingCentre (size, size));
}

void Knob::paint (Graphics& g)
{
    auto r = getLocalBounds().toFloat();
    r.removeFromTop ((float) slider.getBottom() + 4.0f);

    g.setColour (tc::textDim);
    g.setFont (theme::label (12.0f));
    g.drawText (labelText, r.removeFromTop (20.0f), Justification::centred, false);

    g.setColour (tc::text);
    g.setFont (theme::sans (15.0f));
    g.drawText (parameter->getCurrentValueAsText(), r.removeFromTop (24.0f), Justification::centred, false);
}

//==============================================================================
void PresetBar::setPreset (int n, const String& presetName)
{
    number = n;
    name = presetName;
    repaint();
}

Rectangle<float> PresetBar::leftArrow() const { return getLocalBounds().toFloat().removeFromLeft (44.0f); }
Rectangle<float> PresetBar::rightArrow() const { return getLocalBounds().toFloat().removeFromRight (44.0f); }

void PresetBar::paint (Graphics& g)
{
    auto r = getLocalBounds().toFloat().reduced (1.0f);
    g.setColour (Colour (0x99101540));
    g.fillRoundedRectangle (r, 10.0f);
    g.setColour (tc::gold.withAlpha (0.35f));
    g.drawRoundedRectangle (r, 10.0f, 1.0f);

    auto chevron = [&g] (Rectangle<float> a, bool left, bool hot)
    {
        const auto c = a.getCentre();
        Path p;
        const float s = 5.0f;
        p.startNewSubPath (c.x + (left ? s * 0.5f : -s * 0.5f), c.y - s);
        p.lineTo (c.x + (left ? -s * 0.5f : s * 0.5f), c.y);
        p.lineTo (c.x + (left ? s * 0.5f : -s * 0.5f), c.y + s);
        g.setColour (tc::gold.withAlpha (hot ? 1.0f : 0.8f));
        g.strokePath (p, PathStrokeType (1.6f, PathStrokeType::curved, PathStrokeType::rounded));
    };
    chevron (leftArrow(), true, hover == -1);
    chevron (rightArrow(), false, hover == 1);

    // "001  Open Grand", centred as one line
    const auto numberText = String (number).paddedLeft ('0', 3);
    const auto numberFont = theme::sans (16.0f);
    const auto nameFont = theme::sans (17.0f, 0.02f);
    const float gap = 12.0f;
    const float w1 = GlyphArrangement::getStringWidth (numberFont, numberText);
    const float w2 = GlyphArrangement::getStringWidth (nameFont, name);
    float x = r.getCentreX() - (w1 + gap + w2) * 0.5f;

    g.setFont (numberFont);
    g.setColour (tc::textDim);
    g.drawText (numberText, Rectangle<float> (x, r.getY(), w1 + 2, r.getHeight()), Justification::centredLeft, false);
    x += w1 + gap;
    g.setFont (nameFont);
    g.setColour (tc::text);
    g.drawText (name, Rectangle<float> (x, r.getY(), w2 + 2, r.getHeight()), Justification::centredLeft, false);
}

void PresetBar::mouseUp (const MouseEvent& e)
{
    if (! onStep)
        return;
    if (leftArrow().contains (e.position)) onStep (-1);
    else if (rightArrow().contains (e.position)) onStep (1);
}

void PresetBar::mouseMove (const MouseEvent& e)
{
    const int h = leftArrow().contains (e.position) ? -1 : rightArrow().contains (e.position) ? 1 : 0;
    if (h != hover) { hover = h; repaint(); }
    setMouseCursor (h != 0 ? MouseCursor::PointingHandCursor : MouseCursor::NormalCursor);
}

void PresetBar::mouseExit (const MouseEvent&)
{
    hover = 0;
    repaint();
}

//==============================================================================
void StatusLine::setItems (StringPairArray items)
{
    StringArray newKeys = items.getAllKeys(), newValues = items.getAllValues();
    if (newKeys == keys && newValues == values)
        return;
    keys = newKeys;
    values = newValues;
    repaint();
}

void StatusLine::paint (Graphics& g)
{
    // Laid out right to left so the block hugs the right edge
    auto r = getLocalBounds().toFloat();
    const auto keyFont = theme::label (11.5f);
    const auto valueFont = theme::sans (14.0f);

    for (int i = keys.size(); --i >= 0;)
    {
        const float wv = GlyphArrangement::getStringWidth (valueFont, values[i]) + 2;
        const float wk = GlyphArrangement::getStringWidth (keyFont, keys[i]) + 2;
        g.setFont (valueFont);
        g.setColour (tc::text);
        g.drawText (values[i], r.removeFromRight (wv), Justification::centredRight, false);
        r.removeFromRight (9.0f);
        g.setFont (keyFont);
        g.setColour (tc::textDim);
        g.drawText (keys[i], r.removeFromRight (wk), Justification::centredRight, false);
        r.removeFromRight (34.0f);
    }
}

//==============================================================================
NightKeyboard::NightKeyboard (MidiKeyboardState& keyboardState)
    : MidiKeyboardComponent (keyboardState, MidiKeyboardComponent::horizontalKeyboard)
{
    setColour (MidiKeyboardComponent::shadowColourId, Colours::transparentBlack);
    setColour (MidiKeyboardComponent::upDownButtonBackgroundColourId, tc::panel);
    setColour (MidiKeyboardComponent::upDownButtonArrowColourId, tc::gold);
    setScrollButtonsVisible (false);
    setOpaque (false);
}

void NightKeyboard::drawWhiteNote (int note, Graphics& g, Rectangle<float> area, bool isDown, bool isOver, Colour, Colour)
{
    auto key = area.reduced (1.0f, 0.0f).withTrimmedBottom (1.0f);
    ColourGradient ivory (Colour (0xfff7f2e3), key.getTopLeft(), Colour (0xffd9d3c1), key.getBottomLeft(), false);
    g.setGradientFill (ivory);
    g.fillRoundedRectangle (key, 3.0f);

    if (isDown || isOver)
    {
        g.setColour (tc::gold.withAlpha (isDown ? 0.55f : 0.18f));
        g.fillRoundedRectangle (key, 3.0f);
    }

    if (note % 12 == 0) // label the Cs
    {
        g.setColour (Colour (0xff6c6a86));
        g.setFont (theme::sans (10.0f));
        g.drawText (MidiMessage::getMidiNoteName (note, true, true, 4), key.removeFromBottom (16.0f),
                    Justification::centred, false);
    }
}

void NightKeyboard::drawBlackNote (int, Graphics& g, Rectangle<float> area, bool isDown, bool isOver, Colour)
{
    auto key = area.withTrimmedBottom (2.0f);
    ColourGradient ebony (Colour (0xff232a5e), key.getTopLeft(), Colour (0xff0d1130), key.getBottomLeft(), false);
    g.setGradientFill (ebony);
    g.fillRoundedRectangle (key, 2.5f);
    g.setColour (Colours::white.withAlpha (0.10f));
    g.drawRoundedRectangle (key.reduced (0.5f), 2.5f, 0.8f);

    if (isDown || isOver)
    {
        g.setColour (tc::gold.withAlpha (isDown ? 0.6f : 0.2f));
        g.fillRoundedRectangle (key, 2.5f);
    }
}

} // namespace widgets
