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
SkyLookAndFeel::SkyLookAndFeel()
{
    setColour (TextButton::textColourOffId, tc::ink);
    setColour (TextButton::textColourOnId, tc::ink);
    setColour (TooltipWindow::backgroundColourId, Colour (0xf2141834));
    setColour (TooltipWindow::textColourId, tc::ink);
    setColour (TooltipWindow::outlineColourId, tc::glassBorder);
}

void SkyLookAndFeel::drawPearlKnob (Graphics& g, Rectangle<float> bounds, float angle, bool hot, bool showPointer)
{
    const float r = jmin (bounds.getWidth(), bounds.getHeight()) * 0.5f;
    const auto c = bounds.getCentre();
    const auto disc = Rectangle<float> (r * 2, r * 2).withCentre (c);

    // Soft cast shadow, slightly cool
    for (int i = 6; i >= 1; --i)
    {
        g.setColour (tc::shade.withAlpha (0.07f));
        g.fillEllipse (disc.expanded ((float) i * 0.8f).translated (0, (float) i * 0.9f));
    }

    // Pearl body: bright from the upper left, cool blue-grey underneath, warm bounce light
    // from the sun along the lower rim
    ColourGradient body (tc::pearlLight, c.translated (-r * 0.35f, -r * 0.55f),
                         tc::pearlDark, c.translated (r * 0.5f, r * 0.95f), true);
    g.setGradientFill (body);
    g.fillEllipse (disc);

    ColourGradient bounce (tc::sunGlow.withAlpha (0.0f), c, tc::sunGlow.withAlpha (0.28f), c.translated (0, r), false);
    g.setGradientFill (bounce);
    g.fillEllipse (disc);

    // Glossy highlight
    ColourGradient gloss (Colours::white.withAlpha (0.9f), c.translated (-r * 0.3f, -r * 0.5f),
                          Colours::white.withAlpha (0.0f), c.translated (-r * 0.3f, -r * 0.05f), true);
    g.setGradientFill (gloss);
    g.fillEllipse (disc.reduced (r * 0.12f));

    // Rim
    g.setColour (hot ? tc::amber.withAlpha (0.8f) : Colours::white.withAlpha (0.55f));
    g.drawEllipse (disc.reduced (0.5f), hot ? 1.4f : 1.0f);
    g.setColour (tc::shade.withAlpha (0.35f));
    g.drawEllipse (disc.expanded (0.5f), 0.8f);

    if (! showPointer)
        return;

    // Sunlit pointer: an amber bead with a soft halo
    const auto dot = c.getPointOnCircumference (r * 0.74f, angle);
    const float dr = jmax (2.6f, r * 0.08f);
    ColourGradient halo (tc::amber.withAlpha (0.55f), dot, tc::amber.withAlpha (0.0f), dot.translated (dr * 3.5f, 0), true);
    g.setGradientFill (halo);
    g.fillEllipse (Rectangle<float> (dr * 7, dr * 7).withCentre (dot));
    ColourGradient bead (Colour (0xffffe3a8), dot.translated (-dr * 0.3f, -dr * 0.3f), tc::amber, dot.translated (dr, dr), true);
    g.setGradientFill (bead);
    g.fillEllipse (Rectangle<float> (dr * 2, dr * 2).withCentre (dot));
}

void SkyLookAndFeel::drawRotarySlider (Graphics& g, int x, int y, int width, int height, float pos,
                                       float startAngle, float endAngle, Slider& slider)
{
    const auto bounds = Rectangle<float> ((float) x, (float) y, (float) width, (float) height).reduced (6.0f);
    drawPearlKnob (g, bounds, startAngle + pos * (endAngle - startAngle), slider.isMouseOverOrDragging(), true);
}

void SkyLookAndFeel::drawLinearSlider (Graphics& g, int x, int y, int width, int height, float sliderPos,
                                       float, float, Slider::SliderStyle, Slider& slider)
{
    const auto area = Rectangle<float> ((float) x, (float) y, (float) width, (float) height);
    const float trackH = 6.0f;
    const auto track = area.withSizeKeepingCentre (area.getWidth() - 16.0f, trackH);

    // Groove
    g.setColour (tc::shade.withAlpha (0.45f));
    g.fillRoundedRectangle (track, trackH * 0.5f);
    g.setColour (Colours::white.withAlpha (0.25f));
    g.drawRoundedRectangle (track.expanded (0.5f), trackH * 0.5f, 0.8f);

    // Filled part warms from pale gold to amber
    const auto filled = track.withRight (jlimit (track.getX(), track.getRight(), sliderPos));
    ColourGradient fill (Colour (0xff8a4a6a), filled.getTopLeft(), tc::amber, filled.getTopRight(), false);
    g.setGradientFill (fill);
    g.fillRoundedRectangle (filled, trackH * 0.5f);

    // Thumb: a small pearl
    const float thumb = 18.0f;
    drawPearlKnob (g, Rectangle<float> (thumb, thumb).withCentre ({ sliderPos, track.getCentreY() }), 0.0f,
                   slider.isMouseOverOrDragging(), false);
}

void SkyLookAndFeel::drawButtonBackground (Graphics& g, Button& button, const Colour&, bool isHighlighted, bool isDown)
{
    auto r = button.getLocalBounds().toFloat().reduced (glowMargin);
    const bool pill = button.getProperties()[shapeKey].toString() == "pill";
    const float corner = pill ? r.getHeight() * 0.5f : 8.0f;
    const bool on = button.getToggleState();
    const bool enabled = button.isEnabled();

    if (on)
        theme::drawGlow (g, r, corner, tc::amber, glowMargin * 2.0f - 2.0f);

    float opacity = enabled ? 1.0f : 0.55f;
    if (enabled && isHighlighted) opacity = 1.2f;
    if (enabled && isDown) opacity = 1.4f;
    theme::drawGlassPanel (g, r, corner, jmin (1.0f, opacity));
    if (opacity > 1.0f)
    {
        g.setColour (Colours::white.withAlpha ((opacity - 1.0f) * 0.25f));
        g.fillRoundedRectangle (r, corner);
    }

    if (on)
    {
        g.setColour (tc::amber.withAlpha (0.12f));
        g.fillRoundedRectangle (r, corner);
        g.setColour (tc::amber);
        g.drawRoundedRectangle (r, corner, 1.5f);
    }
}

Font SkyLookAndFeel::getTextButtonFont (TextButton& b, int height)
{
    const bool pill = b.getProperties()[shapeKey].toString() == "pill";
    return pill ? theme::label (jmin (12.5f, height * 0.36f)) : theme::sans (jmin (15.0f, height * 0.38f));
}

void SkyLookAndFeel::drawButtonText (Graphics& g, TextButton& button, bool, bool)
{
    auto area = button.getLocalBounds().toFloat().reduced (glowMargin);
    const float alpha = button.isEnabled() ? 1.0f : 0.5f;

    // A/B button: the active letter is amber
    if (button.getProperties().contains (activeSlotKey))
    {
        const bool b = (int) button.getProperties()[activeSlotKey] == 1;
        AttributedString as;
        as.setJustification (Justification::centred);
        const auto f = theme::sans (15.0f, 0.1f);
        as.append ("A", f, b ? tc::inkDim : tc::amber);
        as.append ("  /  ", f, tc::inkDim);
        as.append ("B", f, b ? tc::amber : tc::inkDim);
        as.draw (g, area);
        return;
    }

    const bool pill = button.getProperties()[shapeKey].toString() == "pill";
    const auto text = pill ? button.getButtonText().toUpperCase() : button.getButtonText();
    const auto badge = button.getProperties()[badgeKey].toString();

    g.setColour ((button.getToggleState() && pill ? tc::amber.brighter (0.2f) : tc::ink).withMultipliedAlpha (alpha));
    g.setFont (getTextButtonFont (button, (int) area.getHeight()));

    if (badge.isEmpty())
    {
        g.drawText (text, area, Justification::centred, false);
        return;
    }

    // Label with a small badge underneath, e.g. "soon"
    auto top = area.removeFromTop (area.getHeight() * 0.6f);
    g.drawText (text, top.withTrimmedTop (2.0f), Justification::centredBottom, false);
    g.setColour (tc::amber.withAlpha (0.9f));
    g.setFont (theme::label (8.5f));
    g.drawText (badge.toUpperCase(), area.withTrimmedBottom (2.0f), Justification::centredTop, false);
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
    const int size = jmin (r.getWidth(), r.getHeight() - 44);
    slider.setBounds (r.removeFromTop (size).withSizeKeepingCentre (size, size));
}

void Knob::paint (Graphics& g)
{
    auto r = getLocalBounds().toFloat();
    r.removeFromTop ((float) slider.getBottom() + 2.0f);

    g.setColour (tc::inkDim);
    g.setFont (theme::label (11.5f));
    g.drawText (labelText, r.removeFromTop (18.0f), Justification::centred, false);

    g.setColour (tc::ink);
    g.setFont (theme::sans (14.5f));
    g.drawText (parameter->getCurrentValueAsText(), r.removeFromTop (22.0f), Justification::centred, false);
}

//==============================================================================
KnobSlot::KnobSlot()
{
    setInterceptsMouseClicks (false, false);
}

void KnobSlot::paint (Graphics& g)
{
    // Same footprint as a Knob: a faint pearl and a dashed placeholder for its label
    auto r = getLocalBounds().toFloat();
    const float size = jmin (r.getWidth(), r.getHeight() - 44.0f);
    {
        Graphics::ScopedSaveState s (g);
        g.setOpacity (0.38f);
        SkyLookAndFeel::drawPearlKnob (g, r.removeFromTop (size).withSizeKeepingCentre (size, size).reduced (6.0f), 0, false, false);
    }
    r.removeFromTop (2.0f);
    const auto line = r.removeFromTop (18.0f).withSizeKeepingCentre (34.0f, 1.0f);
    g.setColour (tc::inkDim.withAlpha (0.35f));
    g.fillRect (line);
}

ButtonSlot::ButtonSlot()
{
    setInterceptsMouseClicks (false, false);
}

void ButtonSlot::paint (Graphics& g)
{
    const auto r = getLocalBounds().toFloat().reduced (SkyLookAndFeel::glowMargin);
    Path outline;
    outline.addRoundedRectangle (r, r.getHeight() * 0.5f);
    Path dashed;
    const float dashes[] = { 4.0f, 4.0f };
    PathStrokeType (1.0f).createDashedStroke (dashed, outline, dashes, 2);
    g.setColour (Colours::white.withAlpha (0.04f));
    g.fillPath (outline);
    g.setColour (tc::inkDim.withAlpha (0.4f));
    g.fillPath (dashed);
}

//==============================================================================
HorizontalControl::HorizontalControl (AudioProcessorValueTreeState& state, const String& paramId, const String& text)
    : labelText (text.toUpperCase()),
      parameter (state.getParameter (paramId)),
      attachment (state, paramId, slider)
{
    slider.setSliderStyle (Slider::LinearHorizontal);
    slider.setTextBoxStyle (Slider::NoTextBox, false, 0, 0);
    slider.setDoubleClickReturnValue (true, (double) parameter->convertFrom0to1 (parameter->getDefaultValue()));
    slider.onValueChange = [this] { repaint(); };
    slider.setTooltip (parameter->getName (64) + " (double-click to reset)");
    addAndMakeVisible (slider);
}

void HorizontalControl::resized()
{
    auto r = getLocalBounds();
    r.removeFromLeft (92);
    r.removeFromRight (64);
    slider.setBounds (r);
}

void HorizontalControl::paint (Graphics& g)
{
    auto r = getLocalBounds().toFloat();
    g.setColour (tc::inkDim);
    g.setFont (theme::label (11.5f));
    g.drawText (labelText, r.removeFromLeft (92.0f), Justification::centredLeft, false);
    g.setColour (tc::ink);
    g.setFont (theme::sans (14.5f));
    g.drawText (parameter->getCurrentValueAsText(), r.removeFromRight (60.0f), Justification::centredRight, false);
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
    theme::drawGlassPanel (g, r, 10.0f);
    g.setColour (tc::amber.withAlpha (0.45f));
    g.drawRoundedRectangle (r.reduced (0.5f), 10.0f, 1.0f);

    auto chevron = [&g] (Rectangle<float> a, bool left, bool hot)
    {
        const auto c = a.getCentre();
        Path p;
        const float s = 5.0f;
        p.startNewSubPath (c.x + (left ? s * 0.5f : -s * 0.5f), c.y - s);
        p.lineTo (c.x + (left ? -s * 0.5f : s * 0.5f), c.y);
        p.lineTo (c.x + (left ? s * 0.5f : -s * 0.5f), c.y + s);
        g.setColour (tc::amber.brighter (hot ? 0.25f : 0.0f));
        g.strokePath (p, PathStrokeType (1.7f, PathStrokeType::curved, PathStrokeType::rounded));
    };
    chevron (leftArrow(), true, hover == -1);
    chevron (rightArrow(), false, hover == 1);

    // "001  Open Grand", centred as one line
    const auto numberText = String (number).paddedLeft ('0', 3);
    const auto numberFont = theme::sans (15.5f);
    const auto nameFont = theme::sans (17.0f, 0.02f);
    const float gap = 12.0f;
    const float w1 = GlyphArrangement::getStringWidth (numberFont, numberText);
    const float w2 = GlyphArrangement::getStringWidth (nameFont, name);
    float x = r.getCentreX() - (w1 + gap + w2) * 0.5f;

    g.setFont (numberFont);
    g.setColour (tc::inkDim);
    g.drawText (numberText, Rectangle<float> (x, r.getY(), w1 + 2, r.getHeight()), Justification::centredLeft, false);
    x += w1 + gap;
    g.setFont (nameFont);
    g.setColour (tc::ink);
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
    const auto keyFont = theme::label (10.5f);
    const auto valueFont = theme::sans (13.0f);

    for (int i = keys.size(); --i >= 0;)
    {
        const float wv = GlyphArrangement::getStringWidth (valueFont, values[i]) + 2;
        const float wk = GlyphArrangement::getStringWidth (keyFont, keys[i]) + 2;
        g.setFont (valueFont);
        g.setColour (tc::ink);
        g.drawText (values[i], r.removeFromRight (wv), Justification::centredRight, false);
        r.removeFromRight (8.0f);
        g.setFont (keyFont);
        g.setColour (tc::inkDim);
        g.drawText (keys[i], r.removeFromRight (wk), Justification::centredRight, false);
        r.removeFromRight (26.0f);
    }
}

//==============================================================================
GlassKeyboard::GlassKeyboard (MidiKeyboardState& keyboardState)
    : MidiKeyboardComponent (keyboardState, MidiKeyboardComponent::horizontalKeyboard),
      state (keyboardState)
{
    // Let the sky show through: no background, shadow or separator line
    setColour (MidiKeyboardComponent::whiteNoteColourId, Colours::transparentBlack);
    setColour (MidiKeyboardComponent::shadowColourId, Colours::transparentBlack);
    setColour (MidiKeyboardComponent::keySeparatorLineColourId, Colours::transparentBlack);
    setScrollButtonsVisible (false);
    setOpaque (false);

    for (auto& s : struck)
        s = false;
    state.addListener (&driver);
    driver.startTimerHz (60);
}

GlassKeyboard::~GlassKeyboard()
{
    driver.stopTimer();
    state.removeListener (&driver);
}

void GlassKeyboard::GlowDriver::handleNoteOn (MidiKeyboardState*, int, int note, float)
{
    // May be called from the audio thread: just flag it for the timer
    owner.struck[(size_t) (note & 127)].store (true, std::memory_order_relaxed);
}

void GlassKeyboard::setGlow (int note, float level)
{
    glow[(size_t) (note & 127)] = level;
    repaint();
    if (onGlowChanged)
        onGlowChanged();
}

void GlassKeyboard::updateGlow()
{
    const float decay = std::exp (-1.0f / (60.0f * fadeSeconds));
    bool changed = false;

    for (int note = getRangeStart(); note <= getRangeEnd(); ++note)
    {
        const auto i = (size_t) note;
        const bool held = state.isNoteOnForChannels (0xffff, note);
        const bool hit = struck[i].exchange (false, std::memory_order_relaxed);
        float level = glow[i];

        if (held || hit)
            level = 1.0f;                                  // lit while held (or just tapped)
        else if (wasHeld[i])
            level = jmin (level, remnantLevel);            // released: drop to the remnant
        else if (level > 0.0f)
            level = level * decay < 0.004f ? 0.0f : level * decay; // fade away

        wasHeld[i] = held || hit;
        if (level != glow[i])
        {
            glow[i] = level;
            changed = true;
        }
    }

    if (changed)
    {
        repaint();
        if (onGlowChanged)
            onGlowChanged();
    }
}

void GlassKeyboard::setStyle (KeybedStyle newStyle)
{
    style = newStyle;
    repaint();
}

static void drawCLabel (Graphics& g, int note, Rectangle<float> key, Colour colour)
{
    if (note % 12 != 0)
        return;
    g.setColour (colour);
    g.setFont (theme::sans (10.0f));
    g.drawText (MidiMessage::getMidiNoteName (note, true, true, 4), key.removeFromBottom (18.0f), Justification::centred, false);
}

void GlassKeyboard::drawWhiteNote (int note, Graphics& g, Rectangle<float> area, bool isDown, bool isOver, Colour, Colour)
{
    const float lit = glow[(size_t) note];
    auto key = area.reduced (1.5f, 0.0f).withTrimmedBottom (1.0f);
    const float corner = 5.0f;

    if (style == KeybedStyle::smokedGlass)
    {
        // Smoked glass: the sky shows through, darkened, and each key's front edge catches
        // the sunset glow from below
        ColourGradient glass (Colour (0xff0c1030).withAlpha (0.30f), key.getTopLeft(),
                              Colour (0xff1a1430).withAlpha (0.45f), key.getBottomLeft(), false);
        g.setGradientFill (glass);
        g.fillRoundedRectangle (key, corner);

        ColourGradient rim (tc::sunset.withAlpha (0.0f), { 0, key.getBottom() - key.getHeight() * 0.35f },
                            tc::gold.withAlpha (0.35f), { 0, key.getBottom() }, false);
        g.setGradientFill (rim);
        g.fillRoundedRectangle (key, corner);

        if (lit > 0.0f)
        {
            ColourGradient light (tc::peach.withAlpha (0.35f * lit), key.getTopLeft(),
                                  tc::sunset.withAlpha (0.95f * lit), key.getBottomLeft(), false);
            g.setGradientFill (light);
            g.fillRoundedRectangle (key, corner);
        }
        if (isOver && ! isDown)
        {
            g.setColour (Colours::white.withAlpha (0.08f));
            g.fillRoundedRectangle (key, corner);
        }

        ColourGradient edge (Colours::white.withAlpha (0.30f), key.getTopLeft(),
                             tc::gold.withAlpha (0.75f), key.getBottomLeft(), false);
        g.setGradientFill (edge);
        g.drawRoundedRectangle (key.reduced (0.5f), corner, 1.0f + lit * 0.6f);
        g.setColour (Colours::white.withAlpha (0.12f));
        g.fillRect (key.reduced (0.0f, corner).withX (key.getX() + 2.0f).withWidth (1.2f));
        drawCLabel (g, note, key, tc::ink.withAlpha (0.75f));
        return;
    }

    if (style == KeybedStyle::silhouette)
    {
        // Silhouette: solid black keys like the mountain ridge, with only a hairline of
        // sunset light along the top of each one. Lit keys glow from within like the horizon.
        auto k = area.reduced (0.6f, 0.0f).withTrimmedBottom (1.0f);
        ColourGradient body (Colour (0xff0b0b12), k.getTopLeft(), Colour (0xff040407), k.getBottomLeft(), false);
        g.setGradientFill (body);
        g.fillRoundedRectangle (k, 3.0f);

        if (lit > 0.0f)
        {
            ColourGradient light (Colour (0xff8a3a1c).withAlpha (0.55f * lit), k.getTopLeft(),
                                  tc::sunset.withAlpha (0.95f * lit), { 0, k.getBottom() - k.getHeight() * 0.25f }, false);
            light.addColour (1.0, tc::gold.withAlpha (lit));
            g.setGradientFill (light);
            g.fillRoundedRectangle (k, 3.0f);
        }
        if (isOver && ! isDown)
        {
            g.setColour (Colours::white.withAlpha (0.05f));
            g.fillRoundedRectangle (k, 3.0f);
        }

        // Gap line between keys and a warm lip at the front
        g.setColour (Colour (0xff2a2230));
        g.fillRect (k.withX (k.getRight() - 0.8f).withWidth (0.8f));
        g.setColour (tc::gold.withAlpha (0.55f + 0.45f * lit));
        g.fillRoundedRectangle (k.withTop (k.getBottom() - 2.0f), 1.0f);
        drawCLabel (g, note, k.withTrimmedBottom (2.0f), tc::inkDim.withAlpha (0.8f));
        return;
    }

    // Horizon: each key carries the sky itself, from slate at the back through mauve and
    // peach to a sunset-orange front edge. Lit keys brighten towards gold.
    ColourGradient sky (tc::skyMid.withAlpha (0.55f), key.getTopLeft(), tc::sunset.withAlpha (0.85f), key.getBottomLeft(), false);
    sky.addColour (0.45, tc::skyLow.withAlpha (0.60f));
    sky.addColour (0.78, tc::peach.withAlpha (0.75f));
    g.setGradientFill (sky);
    g.fillRoundedRectangle (key, corner);

    if (lit > 0.0f)
    {
        ColourGradient light (Colour (0xffffe6c0).withAlpha (0.45f * lit), key.getTopLeft(),
                              tc::gold.withAlpha (lit), key.getBottomLeft(), false);
        g.setGradientFill (light);
        g.fillRoundedRectangle (key, corner);
    }
    if (isOver && ! isDown)
    {
        g.setColour (Colours::white.withAlpha (0.12f));
        g.fillRoundedRectangle (key, corner);
    }

    g.setColour (Colours::white.withAlpha (0.35f).interpolatedWith (tc::gold, lit * 0.7f));
    g.drawRoundedRectangle (key.reduced (0.5f), corner, 1.0f + lit * 0.6f);
    drawCLabel (g, note, key, tc::ridge.withAlpha (0.7f));
}

void GlassKeyboard::drawBlackNote (int note, Graphics& g, Rectangle<float> area, bool isDown, bool isOver, Colour)
{
    const float lit = glow[(size_t) note];
    auto key = area.withTrimmedBottom (3.0f);
    const float corner = 4.0f;

    if (style == KeybedStyle::silhouette)
    {
        // Glossy black on matte black: told apart by a faint sheen and a cool rim
        ColourGradient body (Colour (0xff17171f), key.getTopLeft(), Colour (0xff000000), key.getBottomLeft(), false);
        g.setGradientFill (body);
        g.fillRoundedRectangle (key, 2.5f);
        if (lit > 0.0f)
        {
            ColourGradient light (tc::sunset.withAlpha (0.25f * lit), key.getTopLeft(), tc::gold.withAlpha (0.95f * lit), key.getBottomLeft(), false);
            g.setGradientFill (light);
            g.fillRoundedRectangle (key, 2.5f);
        }
        if (isOver && ! isDown)
        {
            g.setColour (Colours::white.withAlpha (0.06f));
            g.fillRoundedRectangle (key, 2.5f);
        }
        g.setColour (tc::inkDim.withAlpha (0.35f));
        g.drawRoundedRectangle (key.reduced (0.5f), 2.5f, 0.8f);
        g.setColour (Colours::white.withAlpha (0.10f));
        g.fillRoundedRectangle (key.reduced (key.getWidth() * 0.25f, 3.0f).withHeight (key.getHeight() * 0.5f), 1.5f);
        return;
    }

    // Glass styles: dark indigo glass, nearly night-coloured
    const auto top = style == KeybedStyle::horizon ? tc::skyHigh : Colour (0xff070914);
    ColourGradient glass (top.withAlpha (0.80f), key.getTopLeft(), tc::skyTop.withAlpha (0.92f), key.getBottomLeft(), false);
    g.setGradientFill (glass);
    g.fillRoundedRectangle (key, corner);

    if (lit > 0.0f)
    {
        ColourGradient light (tc::sunset.withAlpha (0.65f * lit), key.getTopLeft(),
                              tc::gold.withAlpha (0.95f * lit), key.getBottomLeft(), false);
        g.setGradientFill (light);
        g.fillRoundedRectangle (key, corner);
    }

    if (isOver && ! isDown)
    {
        g.setColour (Colours::white.withAlpha (0.10f));
        g.fillRoundedRectangle (key, corner);
    }

    // Glass sheen, and a warm catch-light on the front edge
    g.setColour (Colours::white.withAlpha (0.16f + 0.3f * lit));
    g.drawRoundedRectangle (key.reduced (0.6f), corner, 0.9f);
    g.setColour (Colours::white.withAlpha (0.10f));
    g.fillRoundedRectangle (key.reduced (key.getWidth() * 0.22f, 3.0f).withHeight (key.getHeight() * 0.45f), 2.0f);
    g.setColour (tc::gold.withAlpha (0.45f));
    g.fillRoundedRectangle (key.withTop (key.getBottom() - 1.5f).reduced (2.0f, 0.0f), 0.75f);
}

} // namespace widgets
