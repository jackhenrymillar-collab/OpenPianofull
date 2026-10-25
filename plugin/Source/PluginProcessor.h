/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#pragma once

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_audio_utils/juce_audio_utils.h>

#include "openpiano/PianoInstrument.h"

class OpenPianoAudioProcessor final : public juce::AudioProcessor
{
public:
    OpenPianoAudioProcessor();
    ~OpenPianoAudioProcessor() override;

    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
    void processBlock (juce::AudioBuffer<float>&, juce::MidiBuffer&) override;
    using AudioProcessor::processBlock;

    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override { return true; }

    const juce::String getName() const override { return JucePlugin_Name; }
    bool acceptsMidi() const override { return true; }
    bool producesMidi() const override { return false; }
    bool isMidiEffect() const override { return false; }
    double getTailLengthSeconds() const override { return 0.0; }

    // Factory presets
    int getNumPrograms() override;
    int getCurrentProgram() override { return currentProgram; }
    void setCurrentProgram (int index) override;
    const juce::String getProgramName (int index) override;
    void changeProgramName (int, const juce::String&) override {}

    void getStateInformation (juce::MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

    // A/B comparison: two parameter snapshots. Message thread only.
    int getActiveSlot() const { return activeSlot; }
    void switchToSlot (int slot);

    // For the UI
    juce::AudioProcessorValueTreeState& getState() { return state; }
    juce::MidiKeyboardState& getKeyboardState() { return keyboardState; }
    const openpiano::Instrument& getInstrument() const { return piano; }
    float getCpuLoad() const { return cpuLoad.load (std::memory_order_relaxed); }
    bool isSustainDown() const { return sustainShown.load (std::memory_order_relaxed); }
    double getCurrentSampleRate() const { return currentSampleRate.load (std::memory_order_relaxed); }

private:
    void applyEvent (const juce::MidiMessage& message);
    void updateControls();
    void updateSustain();

    juce::AudioProcessorValueTreeState state;
    juce::MidiKeyboardState keyboardState;

    openpiano::PianoInstrument piano;
    bool prepared = false;

    std::atomic<float>* volumeParam = nullptr;
    std::atomic<float>* dynamicsParam = nullptr;
    std::atomic<float>* hardnessParam = nullptr;
    std::atomic<float>* damperParam = nullptr;
    std::atomic<float>* latchParam = nullptr;
    std::atomic<float>* velocitySensParam = nullptr;

    bool midiSustain = false;
    bool sustainApplied = false;
    std::atomic<bool> sustainShown { false };

    std::atomic<float> cpuLoad { 0.0f };
    std::atomic<double> currentSampleRate { 0.0 };

    int currentProgram = 0;
    int activeSlot = 0;
    std::array<juce::ValueTree, 2> slots;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (OpenPianoAudioProcessor)
};
