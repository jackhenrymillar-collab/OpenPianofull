/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_audio_utils/juce_audio_utils.h>

#include <array>
#include <atomic>

// Dusk-theme controls
namespace widgets
{

class SkyLookAndFeel final : public juce::LookAndFeel_V4
{
public:
    SkyLookAndFeel();

    void drawRotarySlider (juce::Graphics&, int x, int y, int width, int height, float sliderPos,
                           float startAngle, float endAngle, juce::Slider&) override;
    void drawLinearSlider (juce::Graphics&, int x, int y, int width, int height, float sliderPos,
                           float minSliderPos, float maxSliderPos, juce::Slider::SliderStyle, juce::Slider&) override;

    void drawButtonBackground (juce::Graphics&, juce::Button&, const juce::Colour&,
                               bool isHighlighted, bool isDown) override;
    void drawButtonText (juce::Graphics&, juce::TextButton&, bool isHighlighted, bool isDown) override;
    juce::Font getTextButtonFont (juce::TextButton&, int buttonHeight) override;

    // Pearl knob, also used for reserved knob slots (dimmed, no pointer)
    static void drawPearlKnob (juce::Graphics&, juce::Rectangle<float> bounds, float angle, bool hot, bool showPointer);

    // Button shapes, set with button.getProperties().set (shapeKey, ...)
    static constexpr auto shapeKey = "shape";       // "pill" or "tile" (default)
    static constexpr auto badgeKey = "badge";       // small text under the label, e.g. "soon"
    static constexpr auto activeSlotKey = "abSlot"; // marks the A/B button: 0 = A, 1 = B

    // Space kept around every button for the selected glow
    static constexpr float glowMargin = 6.0f;
};

// A pearl knob with a spaced-caps label and its value underneath.
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

// An empty knob position reserved for a future control.
class KnobSlot final : public juce::Component
{
public:
    KnobSlot();
    void paint (juce::Graphics&) override;
};

// An empty button position reserved for a future control.
class ButtonSlot final : public juce::Component
{
public:
    ButtonSlot();
    void paint (juce::Graphics&) override;
};

// A horizontal slider with its label on the left and value on the right.
class HorizontalControl final : public juce::Component
{
public:
    HorizontalControl (juce::AudioProcessorValueTreeState& state, const juce::String& paramId, const juce::String& labelText);

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

// "KEY value   KEY value", right aligned
class StatusLine final : public juce::Component
{
public:
    void setItems (juce::StringPairArray items);
    void paint (juce::Graphics&) override;

private:
    juce::StringArray keys, values;
};

// The three keybed looks offered for the dusk theme
enum class KeybedStyle { smokedGlass, silhouette, horizon };

// A see-through glass keyboard. Played keys light up with sunlight and, once released,
// leave a glowing remnant that fades away. Notes from the host light it too.
class GlassKeyboard final : public juce::MidiKeyboardComponent
{
public:
    explicit GlassKeyboard (juce::MidiKeyboardState& keyboardState);
    ~GlassKeyboard() override;

    // 0..1: how brightly a key is lit
    float glowLevel (int note) const { return glow[(size_t) (note & 127)]; }

    // For previews: light a key as if it had just been released
    void setGlow (int note, float level);

    void setStyle (KeybedStyle);
    KeybedStyle getStyle() const { return style; }

    // Called whenever the glow changes, so the parent can repaint its bloom
    std::function<void()> onGlowChanged;

    void drawWhiteNote (int midiNoteNumber, juce::Graphics&, juce::Rectangle<float> area,
                        bool isDown, bool isOver, juce::Colour lineColour, juce::Colour textColour) override;
    void drawBlackNote (int midiNoteNumber, juce::Graphics&, juce::Rectangle<float> area,
                        bool isDown, bool isOver, juce::Colour noteFillColour) override;

    // Time constant of the fading remnant, in seconds
    static constexpr float fadeSeconds = 1.3f;
    // Brightness the remnant starts at when a key is released
    static constexpr float remnantLevel = 0.75f;

private:
    // MidiKeyboardComponent is already a keyboard listener and a timer, so the glow
    // animation gets its own small helper for both.
    struct GlowDriver final : juce::MidiKeyboardState::Listener, juce::Timer
    {
        explicit GlowDriver (GlassKeyboard& k) : owner (k) {}
        void handleNoteOn (juce::MidiKeyboardState*, int, int note, float) override;
        void handleNoteOff (juce::MidiKeyboardState*, int, int, float) override {}
        void timerCallback() override { owner.updateGlow(); }
        GlassKeyboard& owner;
    };

    void updateGlow();

    juce::MidiKeyboardState& state;
    KeybedStyle style = KeybedStyle::smokedGlass;
    GlowDriver driver { *this };
    std::array<std::atomic<bool>, 128> struck {};
    std::array<float, 128> glow {};
    std::array<bool, 128> wasHeld {};
};

} // namespace widgets
