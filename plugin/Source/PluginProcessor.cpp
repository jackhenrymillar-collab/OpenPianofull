/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"
#include "Parameters.h"

OpenPianoAudioProcessor::OpenPianoAudioProcessor()
    : AudioProcessor (BusesProperties().withOutput ("Output", juce::AudioChannelSet::stereo(), true)),
      state (*this, nullptr, "OpenPiano", params::createLayout())
{
    volumeParam = state.getRawParameterValue (params::volume);
    dynamicsParam = state.getRawParameterValue (params::dynamics);
    hardnessParam = state.getRawParameterValue (params::hardness);
    damperParam = state.getRawParameterValue (params::damper);
    actionParam = state.getRawParameterValue (params::actionNoise);
    ageParam = state.getRawParameterValue (params::age);
    roomParam = state.getRawParameterValue (params::room);
    latchParam = state.getRawParameterValue (params::sustainLatch);
    velocitySensParam = state.getRawParameterValue (params::velocitySens);
}

OpenPianoAudioProcessor::~OpenPianoAudioProcessor() = default;

void OpenPianoAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    piano.prepare (sampleRate, samplesPerBlock);
    keyboardState.reset();
    midiSustain = false;
    sustainApplied = false;
    currentSampleRate = sampleRate;
    prepared = true;
}

void OpenPianoAudioProcessor::releaseResources()
{
    piano.reset();
}

bool OpenPianoAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
    const auto out = layouts.getMainOutputChannelSet();
    return out == juce::AudioChannelSet::mono() || out == juce::AudioChannelSet::stereo();
}

void OpenPianoAudioProcessor::updateControls()
{
    openpiano::PianoControls c;
    c.outputGain = juce::Decibels::decibelsToGain (volumeParam->load(), -100.0f);
    c.dynamicsCurve = dynamicsParam->load();
    c.hammerHardness = hardnessParam->load();
    c.damperStrength = damperParam->load();
    c.actionNoise = actionParam->load();
    c.age = ageParam->load();
    c.room = roomParam->load();

    // Snap the defaults exactly, so the default sound matches the regression-tested engine
    if (std::abs (c.dynamicsCurve - 1.0f) < 0.005f) c.dynamicsCurve = 1.0f;
    if (std::abs (c.hammerHardness - 1.0f) < 0.005f) c.hammerHardness = 1.0f;
    if (std::abs (c.damperStrength - 1.0f) < 0.005f) c.damperStrength = 1.0f;
    if (std::abs (volumeParam->load()) < 0.05f) c.outputGain = 1.0f;

    piano.setControls (c);
}

void OpenPianoAudioProcessor::updateSustain()
{
    const bool down = midiSustain || latchParam->load() > 0.5f;
    if (down != sustainApplied)
    {
        piano.setSustain (down);
        sustainApplied = down;
    }
    sustainShown.store (down, std::memory_order_relaxed);
}

void OpenPianoAudioProcessor::applyEvent (const juce::MidiMessage& m)
{
    if (m.isNoteOn())
    {
        // Velocity sensitivity pulls every velocity towards a fixed medium-loud touch
        int velocity = m.getVelocity();
        const float sens = velocitySensParam->load();
        if (sens < 0.995f)
            velocity = juce::jlimit (1, 127, juce::roundToInt (96.0f + (velocity - 96.0f) * sens));
        piano.noteOn (m.getNoteNumber(), velocity);
    }
    else if (m.isNoteOff())
        piano.noteOff (m.getNoteNumber());
    else if (m.isSustainPedalOn())
        midiSustain = true;
    else if (m.isSustainPedalOff())
        midiSustain = false;
    else if (m.isSostenutoPedalOn())
        piano.setSostenuto (true);
    else if (m.isSostenutoPedalOff())
        piano.setSostenuto (false);
    else if (m.isController() && m.getControllerNumber() == 67) // soft pedal (una corda), half pedalling allowed
        piano.setSoftPedal (m.getControllerValue() / 127.0f);
    else if (m.isAllNotesOff() || m.isAllSoundOff() || m.isResetAllControllers())
    {
        if (m.isAllSoundOff())
            piano.reset();
        for (int n = 0; n < 128; ++n)
            piano.noteOff (n);
        midiSustain = false;
        if (m.isResetAllControllers())
        {
            piano.setSostenuto (false);
            piano.setSoftPedal (0.0f);
        }
    }
    updateSustain();
}

void OpenPianoAudioProcessor::processBlock (juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midi)
{
    juce::ScopedNoDenormals noDenormals;
    const auto startTicks = juce::Time::getHighResolutionTicks();
    const int numSamples = buffer.getNumSamples();
    const int numChannels = buffer.getNumChannels();

    if (! prepared || numChannels == 0)
    {
        buffer.clear();
        return;
    }

    // Add notes played on the on-screen keyboard
    keyboardState.processNextMidiBuffer (midi, 0, numSamples, true);

    updateControls();
    updateSustain();

    // Render up to each MIDI event, then apply it (sample-accurate timing)
    std::array<float*, 2> channels {};
    int pos = 0;
    auto renderTo = [&] (int end)
    {
        if (end <= pos)
            return;
        for (int ch = 0; ch < std::min (numChannels, 2); ++ch)
            channels[(size_t) ch] = buffer.getWritePointer (ch, pos);
        piano.process (channels.data(), std::min (numChannels, 2), end - pos);
        pos = end;
    };

    for (const auto metadata : midi)
    {
        renderTo (juce::jlimit (0, numSamples, metadata.samplePosition));
        applyEvent (metadata.getMessage());
    }
    renderTo (numSamples);

    for (int ch = 2; ch < numChannels; ++ch)
        buffer.clear (ch, 0, numSamples);

    // CPU load: time spent here relative to the block's duration, smoothed
    const double seconds = juce::Time::highResolutionTicksToSeconds (juce::Time::getHighResolutionTicks() - startTicks);
    const double blockSeconds = numSamples / juce::jmax (1.0, getSampleRate());
    const float load = (float) (seconds / blockSeconds);
    cpuLoad.store (cpuLoad.load (std::memory_order_relaxed) * 0.9f + load * 0.1f, std::memory_order_relaxed);
}

//==============================================================================
int OpenPianoAudioProcessor::getNumPrograms()
{
    return (int) params::factoryPresets().size();
}

const juce::String OpenPianoAudioProcessor::getProgramName (int index)
{
    const auto& presets = params::factoryPresets();
    return juce::isPositiveAndBelow (index, (int) presets.size()) ? juce::String (presets[(size_t) index].name) : juce::String();
}

void OpenPianoAudioProcessor::setCurrentProgram (int index)
{
    const auto& presets = params::factoryPresets();
    if (! juce::isPositiveAndBelow (index, (int) presets.size()))
        return;

    currentProgram = index;
    const auto& p = presets[(size_t) index];
    auto set = [this] (const char* id, float value)
    {
        if (auto* param = state.getParameter (id))
            param->setValueNotifyingHost (param->convertTo0to1 (value));
    };
    set (params::volume, p.volumeDb);
    set (params::dynamics, p.dynamics);
    set (params::hardness, p.hardness);
    set (params::damper, p.damper);
}

void OpenPianoAudioProcessor::switchToSlot (int slot)
{
    if (slot == activeSlot || ! juce::isPositiveAndBelow (slot, 2))
        return;

    auto current = state.copyState();
    current.setProperty ("program", currentProgram, nullptr);
    slots[(size_t) activeSlot] = current;

    // First switch to an empty slot starts it as a copy of the current sound
    if (slots[(size_t) slot].isValid())
    {
        state.replaceState (slots[(size_t) slot]);
        currentProgram = slots[(size_t) slot].getProperty ("program", currentProgram);
    }

    activeSlot = slot;
}

void OpenPianoAudioProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    auto tree = state.copyState();
    tree.setProperty ("program", currentProgram, nullptr);
    if (auto xml = tree.createXml())
        copyXmlToBinary (*xml, destData);
}

void OpenPianoAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    if (auto xml = getXmlFromBinary (data, sizeInBytes))
    {
        if (xml->hasTagName (state.state.getType()))
        {
            auto tree = juce::ValueTree::fromXml (*xml);
            currentProgram = tree.getProperty ("program", 0);
            state.replaceState (tree);
        }
    }
}

juce::AudioProcessorEditor* OpenPianoAudioProcessor::createEditor()
{
    return new OpenPianoAudioProcessorEditor (*this);
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new OpenPianoAudioProcessor();
}
