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

    static constexpr int kWidth = 1040;
    static constexpr int kHeight = 720;

private:
    void timerCallback() override;
    void refreshPreset();
    void refreshAB();

    OpenPianoAudioProcessor& processor;
    widgets::NightLookAndFeel lookAndFeel;
    juce::TooltipWindow tooltips { this, 600 };

    juce::Image sky;

    widgets::PresetBar presetBar;
    juce::TextButton abButton { "A / B" };

    // Instrument modes. Only Piano exists so far.
    juce::TextButton pianoMode { "Piano" }, stringsMode { "Strings" }, rhodesMode { "Rhodes" };

    widgets::Knob volume, dynamics, hardness, damper;

    juce::TextButton sustainButton { "Sustain" };
    juce::AudioProcessorValueTreeState::ButtonAttachment sustainAttachment;
    widgets::StatusLine status;

    widgets::NightKeyboard keyboard;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (OpenPianoAudioProcessorEditor)
};
