/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_audio_utils/juce_audio_utils.h>

// Overcast-style controls
namespace widgets
{

class NightLookAndFeel final : public juce::LookAndFeel_V4
{
public:
    NightLookAndFeel();

    void drawRotarySlider (juce::Graphics&, int x, int y, int width, int height, float sliderPos,
                           float startAngle, float endAngle, juce::Slider&) override;

    void drawButtonBackground (juce::Graphics&, juce::Button&, const juce::Colour&,
                               bool isHighlighted, bool isDown) override;
    void drawButtonText (juce::Graphics&, juce::TextButton&, bool isHighlighted, bool isDown) override;

    juce::Font getTextButtonFont (juce::TextButton&, int buttonHeight) override;

    // Button shapes, set with button.getProperties().set (shapeKey, ...)
    static constexpr auto shapeKey = "shape";   // "pill" or "tile" (default)
    static constexpr auto badgeKey = "badge";   // small text under the label, e.g. "soon"
    static constexpr auto activeSlotKey = "abSlot"; // marks the A/B button: 0 = A, 1 = B

    // Space kept around every button for the selected glow
    static constexpr float glowMargin = 6.0f;
};

// A big round knob with a spaced-caps label and its value underneath.
class Knob final : public juce::Component
{
public:
    Knob (juce::AudioProcessorValueTreeState& state, const juce::String& paramId, const juce::String& labelText);

    void resized() override;
    void paint (juce::Graphics&) override;

private:
    juce::Slider slider;
    juce::String labelText;
    juce::RangedAudioParameter* parameter;
    juce::AudioProcessorValueTreeState::SliderAttachment attachment;
};

// "‹  001  Open Grand  ›"
class PresetBar final : public juce::Component
{
public:
    std::function<void (int delta)> onStep;
    void setPreset (int number, const juce::String& name);

    void paint (juce::Graphics&) override;
    void mouseUp (const juce::MouseEvent&) override;
    void mouseMove (const juce::MouseEvent&) override;
    void mouseExit (const juce::MouseEvent&) override;

private:
    juce::Rectangle<float> leftArrow() const;
    juce::Rectangle<float> rightArrow() const;

    int number = 1;
    juce::String name;
    int hover = 0; // -1 left, 1 right
};

// "KEY value   KEY value" in the bottom right
class StatusLine final : public juce::Component
{
public:
    void setItems (juce::StringPairArray items);
    void paint (juce::Graphics&) override;

private:
    juce::StringArray keys, values;
};

// The on-screen keyboard, drawn as ivory and ebony keys with a gold glow on pressed keys.
class NightKeyboard final : public juce::MidiKeyboardComponent
{
public:
    NightKeyboard (juce::MidiKeyboardState& keyboardState);

    void drawWhiteNote (int midiNoteNumber, juce::Graphics&, juce::Rectangle<float> area,
                        bool isDown, bool isOver, juce::Colour lineColour, juce::Colour textColour) override;
    void drawBlackNote (int midiNoteNumber, juce::Graphics&, juce::Rectangle<float> area,
                        bool isDown, bool isOver, juce::Colour noteFillColour) override;
};

} // namespace widgets
