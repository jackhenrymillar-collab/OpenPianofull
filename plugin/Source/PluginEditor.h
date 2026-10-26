/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include "PluginProcessor.h"
#include "ui/Widgets.h"

class OpenPianoAudioProcessorEditor final : public juce::AudioProcessorEditor,
                                            private juce::Timer
{
public:
    explicit OpenPianoAudioProcessorEditor (OpenPianoAudioProcessor&);
    ~OpenPianoAudioProcessorEditor() override;

    void paint (juce::Graphics&) override;
    void resized() override;

    // For screenshots: light some keys as if just played
    void previewGlow();

    static constexpr int kWidth = 1040;
    static constexpr int kHeight = 720;

private:
    void timerCallback() override;
    void refreshPreset();
    void refreshAB();
    juce::Rectangle<float> deckBounds() const;

    OpenPianoAudioProcessor& processor;
    widgets::SkyLookAndFeel lookAndFeel;
    juce::TooltipWindow tooltips { this, 600 };

    juce::Image skyImage;

    // Top: instrument modes and the keyboard. Only Piano exists so far.
    juce::TextButton pianoMode { "Piano" }, stringsMode { "Strings" }, rhodesMode { "Rhodes" };
    widgets::GlassKeyboard keyboard;

    // Bottom deck
    widgets::PresetBar presetBar;
    juce::TextButton abButton { "A / B" };
    widgets::ButtonSlot saveSlot;

    widgets::Knob volume, dynamics, hardness, damper;
    widgets::Knob action; // first of the smaller knobs
    widgets::Knob age;    // second
    widgets::Knob room;   // third
    std::array<widgets::KnobSlot, 4> smallKnobSlots;

    widgets::HorizontalControl velocity;
    juce::TextButton sustainButton { "Sustain" };
    juce::AudioProcessorValueTreeState::ButtonAttachment sustainAttachment;
    std::array<widgets::ButtonSlot, 3> buttonSlots;

    widgets::StatusLine status;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (OpenPianoAudioProcessorEditor)
};
