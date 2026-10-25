/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#include "PluginEditor.h"
#include "Parameters.h"
#include "ui/NightSky.h"
#include "ui/Theme.h"

using namespace juce;
namespace tc = theme::colours;

namespace
{
    const Point<float> moonCentre { 887.0f, 139.0f };
    constexpr float moonRadius = 74.0f;
}

OpenPianoAudioProcessorEditor::OpenPianoAudioProcessorEditor (OpenPianoAudioProcessor& p)
    : AudioProcessorEditor (&p),
      processor (p),
      volume (p.getState(), params::volume, "Volume"),
      dynamics (p.getState(), params::dynamics, "Dynamics"),
      hardness (p.getState(), params::hardness, "Hammer"),
      damper (p.getState(), params::damper, "Damper"),
      sustainAttachment (p.getState(), params::sustainLatch, sustainButton),
      keyboard (p.getKeyboardState())
{
    setLookAndFeel (&lookAndFeel);

    presetBar.onStep = [this] (int delta)
    {
        const int n = processor.getNumPrograms();
        processor.setCurrentProgram ((processor.getCurrentProgram() + delta + n) % n);
        refreshPreset();
    };
    addAndMakeVisible (presetBar);

    abButton.getProperties().set (widgets::NightLookAndFeel::shapeKey, "tile");
    abButton.setTooltip ("Compare two settings. Click to switch between A and B.");
    abButton.onClick = [this]
    {
        processor.switchToSlot (1 - processor.getActiveSlot());
        refreshAB();
        refreshPreset();
    };
    addAndMakeVisible (abButton);

    for (auto* b : { &pianoMode, &stringsMode, &rhodesMode })
        addAndMakeVisible (b);
    pianoMode.setToggleState (true, dontSendNotification);
    pianoMode.setTooltip ("Hammered strings, simulated in real time");
    for (auto* b : { &stringsMode, &rhodesMode })
    {
        b->setEnabled (false);
        b->getProperties().set (widgets::NightLookAndFeel::badgeKey, "soon");
    }

    for (auto* k : { &volume, &dynamics, &hardness, &damper })
        addAndMakeVisible (k);

    sustainButton.setClickingTogglesState (true);
    sustainButton.getProperties().set (widgets::NightLookAndFeel::shapeKey, "pill");
    sustainButton.setTooltip ("Latch the sustain pedal. It also lights up when your pedal is down.");
    addAndMakeVisible (sustainButton);

    addAndMakeVisible (status);

    const auto& instrument = processor.getInstrument();
    keyboard.setAvailableRange (instrument.lowestNote(), instrument.highestNote());
    keyboard.setOctaveForMiddleC (4);
    keyboard.setKeyPressBaseOctave (4);
    addAndMakeVisible (keyboard);

    setSize (kWidth, kHeight);
    refreshPreset();
    refreshAB();
    timerCallback();
    startTimerHz (15);
}

OpenPianoAudioProcessorEditor::~OpenPianoAudioProcessorEditor()
{
    setLookAndFeel (nullptr);
}

void OpenPianoAudioProcessorEditor::refreshPreset()
{
    const int index = processor.getCurrentProgram();
    presetBar.setPreset (index + 1, processor.getProgramName (index));
}

void OpenPianoAudioProcessorEditor::refreshAB()
{
    abButton.getProperties().set (widgets::NightLookAndFeel::activeSlotKey, processor.getActiveSlot());
    abButton.repaint();
    repaint (abButton.getBounds().expanded (12));
}

void OpenPianoAudioProcessorEditor::timerCallback()
{
    // The latch button shows the effective pedal state (MIDI pedal or latch)
    const bool latched = processor.getState().getRawParameterValue (params::sustainLatch)->load() > 0.5f;
    const bool down = processor.isSustainDown();
    if (sustainButton.getToggleState() != (latched || down))
        sustainButton.setToggleState (latched || down, dontSendNotification);

    const auto& instrument = processor.getInstrument();
    StringPairArray items;
    items.set ("RANGE", MidiMessage::getMidiNoteName (instrument.lowestNote(), true, true, 4) + String::fromUTF8 ("\xe2\x80\x93")
                            + MidiMessage::getMidiNoteName (instrument.highestNote(), true, true, 4));
    items.set ("STRINGS", String (instrument.activeVoiceCount()));
    items.set ("CPU", String (roundToInt (processor.getCpuLoad() * 100.0f)) + "%");
    const double sr = processor.getCurrentSampleRate();
    items.set ("RATE", sr > 0 ? String (sr / 1000.0, 1) + " kHz" : String ("-"));
    status.setItems (items);
}

void OpenPianoAudioProcessorEditor::paint (Graphics& g)
{
    const float scale = (float) g.getInternalContext().getPhysicalPixelScaleFactor();
    if (sky.isNull() || sky.getWidth() != roundToInt (getWidth() * scale))
    {
        sky = Image (Image::RGB, roundToInt (getWidth() * scale), roundToInt (getHeight() * scale), true);
        Graphics sg (sky);
        sg.addTransform (AffineTransform::scale (scale));
        nightsky::paint (sg, getLocalBounds().toFloat(), { moonCentre, moonRadius });
    }
    g.drawImage (sky, getLocalBounds().toFloat());

    // Wordmark
    nightsky::paintCrescent (g, { 44.0f, 34.0f, 26.0f, 26.0f }, tc::moon);
    g.setColour (tc::text);
    g.setFont (theme::serif (34.0f, 0.32f));
    g.drawText (String (JucePlugin_Name).toUpperCase(), Rectangle<float> (80.0f, 26.0f, 320.0f, 42.0f), Justification::centredLeft, false);

    // Mode heading
    g.setColour (tc::textDim);
    g.setFont (theme::label (12.0f));
    g.drawText ("INSTRUMENT", Rectangle<float> (42.0f, 126.0f, 300.0f, 18.0f), Justification::centredLeft, false);

    g.setColour (tc::text);
    g.setFont (theme::serif (56.0f));
    g.drawText ("Piano", Rectangle<float> (40.0f, 146.0f, 500.0f, 62.0f), Justification::centredLeft, false);

    g.setColour (tc::text.withAlpha (0.85f));
    g.setFont (theme::sans (16.0f));
    g.drawText ("Hammered stiff strings, solved in real time with finite differences",
                Rectangle<float> (42.0f, 210.0f, 700.0f, 24.0f), Justification::centredLeft, false);

    // Soft panel behind the keyboard
    auto kb = keyboard.getBounds().toFloat().expanded (10.0f, 10.0f);
    g.setColour (Colour (0x700b0f2c));
    g.fillRoundedRectangle (kb, 12.0f);
    g.setColour (tc::panelBorder.withAlpha (0.12f));
    g.drawRoundedRectangle (kb, 12.0f, 1.0f);
}

void OpenPianoAudioProcessorEditor::resized()
{
    presetBar.setBounds (404, 26, 356, 44);
    const int m = (int) widgets::NightLookAndFeel::glowMargin;
    abButton.setBounds (Rectangle<int> (kWidth - 40 - 84, 31, 84, 34).expanded (m));

    const int modeW = 112, modeH = 46, modeGap = 10;
    int x = 40;
    for (auto* b : { &pianoMode, &stringsMode, &rhodesMode })
    {
        b->setBounds (Rectangle<int> (x, 258, modeW, modeH).expanded (m));
        x += modeW + modeGap;
    }

    // Four knobs, evenly spread
    const int knobW = 150, knobH = 150;
    const int spacing = (kWidth - 80) / 4;
    int i = 0;
    for (auto* k : { &volume, &dynamics, &hardness, &damper })
    {
        const int cx = 40 + spacing * i + spacing / 2;
        k->setBounds (cx - knobW / 2, 330, knobW, knobH);
        ++i;
    }

    sustainButton.setBounds (Rectangle<int> (40, 512, 130, 40).expanded (m));
    status.setBounds (380, 512, kWidth - 40 - 380, 40);

    keyboard.setBounds (40, 584, kWidth - 80, 112);
    const int whiteKeys = [this]
    {
        int n = 0;
        for (int note = keyboard.getRangeStart(); note <= keyboard.getRangeEnd(); ++note)
            n += MidiMessage::isMidiNoteBlack (note) ? 0 : 1;
        return n;
    }();
    keyboard.setKeyWidth ((float) keyboard.getWidth() / (float) jmax (1, whiteKeys));
}
