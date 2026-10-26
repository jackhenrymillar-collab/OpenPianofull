/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

#include "PluginEditor.h"
#include "Parameters.h"
#include "ui/Sky.h"
#include "ui/Theme.h"

using namespace juce;
namespace tc = theme::colours;
using LAF = widgets::SkyLookAndFeel;

OpenPianoAudioProcessorEditor::OpenPianoAudioProcessorEditor (OpenPianoAudioProcessor& p)
    : AudioProcessorEditor (&p),
      processor (p),
      keyboard (p.getKeyboardState()),
      volume (p.getState(), params::volume, "Volume"),
      dynamics (p.getState(), params::dynamics, "Dynamics"),
      hardness (p.getState(), params::hardness, "Hammer"),
      damper (p.getState(), params::damper, "Damper"),
      action (p.getState(), params::actionNoise, "Action"),
      age (p.getState(), params::age, "Age"),
      velocity (p.getState(), params::velocitySens, "Velocity"),
      sustainAttachment (p.getState(), params::sustainLatch, sustainButton)
{
    setLookAndFeel (&lookAndFeel);

    // Modes
    for (auto* b : { &pianoMode, &stringsMode, &rhodesMode })
        addAndMakeVisible (b);
    pianoMode.setToggleState (true, dontSendNotification);
    pianoMode.setTooltip ("Hammered strings, simulated in real time");
    for (auto* b : { &stringsMode, &rhodesMode })
    {
        b->setEnabled (false);
        b->getProperties().set (LAF::badgeKey, "soon");
    }

    // Keyboard
    const auto& instrument = processor.getInstrument();
    keyboard.setAvailableRange (instrument.lowestNote(), instrument.highestNote());
    keyboard.setOctaveForMiddleC (4);
    keyboard.setKeyPressBaseOctave (4);
    keyboard.onGlowChanged = [this] { repaint (keyboard.getBounds().expanded (70, 50)); };
    addAndMakeVisible (keyboard);

    // Presets and A/B
    presetBar.onStep = [this] (int delta)
    {
        const int n = processor.getNumPrograms();
        processor.setCurrentProgram ((processor.getCurrentProgram() + delta + n) % n);
        refreshPreset();
    };
    addAndMakeVisible (presetBar);

    abButton.setTooltip ("Compare two settings. Click to switch between A and B.");
    abButton.onClick = [this]
    {
        processor.switchToSlot (1 - processor.getActiveSlot());
        refreshAB();
        refreshPreset();
    };
    addAndMakeVisible (abButton);
    addAndMakeVisible (saveSlot);

    // Knobs: four main ones, and four smaller positions reserved for later
    for (auto* k : { &volume, &dynamics, &hardness, &damper })
        addAndMakeVisible (k);
    for (auto& s : smallKnobSlots)
        addAndMakeVisible (s);
    smallKnobSlots[0].setVisible (false);
    smallKnobSlots[1].setVisible (false);
    addAndMakeVisible (action);
    addAndMakeVisible (age);

    // Bottom row
    addAndMakeVisible (velocity);

    sustainButton.setClickingTogglesState (true);
    sustainButton.getProperties().set (LAF::shapeKey, "pill");
    sustainButton.setTooltip ("Latch the sustain pedal. It also lights up when your pedal is down.");
    addAndMakeVisible (sustainButton);
    for (auto& s : buttonSlots)
        addAndMakeVisible (s);

    addAndMakeVisible (status);

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

void OpenPianoAudioProcessorEditor::previewGlow()
{
    // A chord being held (C4 E4 G4) and the remnants of an earlier phrase fading out
    for (int n : { 60, 64, 67 })
        keyboard.setGlow (n, 1.0f);
    const std::pair<int, float> fading[] = { { 48, 0.55f }, { 52, 0.38f }, { 55, 0.25f }, { 57, 0.15f }, { 43, 0.08f }, { 70, 0.45f } };
    for (auto [n, level] : fading)
        keyboard.setGlow (n, level);
}

void OpenPianoAudioProcessorEditor::refreshPreset()
{
    const int index = processor.getCurrentProgram();
    presetBar.setPreset (index + 1, processor.getProgramName (index));
}

void OpenPianoAudioProcessorEditor::refreshAB()
{
    abButton.getProperties().set (LAF::activeSlotKey, processor.getActiveSlot());
    abButton.repaint();
}

void OpenPianoAudioProcessorEditor::timerCallback()
{
    // The sustain button shows the effective pedal state (MIDI pedal or latch)
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

Rectangle<float> OpenPianoAudioProcessorEditor::deckBounds() const
{
    return { 24.0f, 318.0f, (float) kWidth - 48.0f, 384.0f };
}

void OpenPianoAudioProcessorEditor::paint (Graphics& g)
{
    const float scale = (float) g.getInternalContext().getPhysicalPixelScaleFactor();
    if (skyImage.isNull() || skyImage.getWidth() != roundToInt (getWidth() * scale))
    {
        skyImage = Image (Image::RGB, roundToInt (getWidth() * scale), roundToInt (getHeight() * scale), true);
        Graphics sg (skyImage);
        sg.addTransform (AffineTransform::scale (scale));
        sky::paint (sg, getLocalBounds().toFloat());
    }
    g.drawImage (skyImage, getLocalBounds().toFloat());

    // Wordmark
    sky::paintSunIcon (g, { 40.0f, 25.0f, 26.0f, 26.0f }, Colours::white);
    g.setColour (Colours::white);
    g.setFont (theme::serif (30.0f, 0.32f));
    g.drawText (String (JucePlugin_Name).toUpperCase(), Rectangle<float> (78.0f, 17.0f, 360.0f, 42.0f), Justification::centredLeft, false);

    // Light spilling from lit keys onto the sky around the keyboard
    {
        const auto origin = keyboard.getPosition().toFloat();
        for (int note = keyboard.getRangeStart(); note <= keyboard.getRangeEnd(); ++note)
        {
            const float lit = keyboard.glowLevel (note);
            if (lit <= 0.0f)
                continue;
            const auto key = keyboard.getRectangleForKey (note) + origin;
            const auto c = Point<float> (key.getCentreX(), key.getBottom() - key.getHeight() * 0.25f);
            const float rx = keyboard.getKeyWidth() * 2.6f;
            Graphics::ScopedSaveState s (g);
            g.addTransform (AffineTransform::scale (1.0f, 1.6f, c.x, c.y));
            ColourGradient bloom (tc::sunGlow.withAlpha (0.55f * lit), c, tc::sunGlow.withAlpha (0.0f), c.translated (rx, 0), true);
            g.setGradientFill (bloom);
            g.fillEllipse (c.x - rx, c.y - rx, rx * 2, rx * 2);
        }
    }

    // A thin glass shelf under the keys
    auto kb = keyboard.getBounds().toFloat();
    auto shelf = Rectangle<float> (kb.getX() - 14.0f, kb.getBottom() + 4.0f, kb.getWidth() + 28.0f, 10.0f);
    theme::drawGlassPanel (g, shelf, 5.0f, 0.9f);

    // The control deck: one large frosted panel
    theme::drawGlassPanel (g, deckBounds(), 18.0f, 0.95f);

    // Divider between the main knobs and the reserved smaller ones
    const auto first = smallKnobSlots.front().getBounds().toFloat();
    g.setColour (Colours::white.withAlpha (0.7f));
    g.fillRect (Rectangle<float> (first.getX() - 22.0f, first.getY() + 6.0f, 1.0f, first.getHeight() - 20.0f));
}

void OpenPianoAudioProcessorEditor::resized()
{
    const int m = (int) LAF::glowMargin;

    // Header: modes on the right
    int x = kWidth - 40 - 3 * 100 - 2 * 8;
    for (auto* b : { &pianoMode, &stringsMode, &rhodesMode })
    {
        b->setBounds (Rectangle<int> (x, 18, 100, 40).expanded (m));
        x += 108;
    }

    // Keyboard centred near the top
    const int keyboardWidth = 920;
    keyboard.setBounds ((kWidth - keyboardWidth) / 2, 86, keyboardWidth, 196);
    int whiteKeys = 0;
    for (int note = keyboard.getRangeStart(); note <= keyboard.getRangeEnd(); ++note)
        whiteKeys += MidiMessage::isMidiNoteBlack (note) ? 0 : 1;
    keyboard.setKeyWidth ((float) keyboard.getWidth() / (float) jmax (1, whiteKeys));

    // Deck, row 1: presets centred, A/B on the right, a reserved button on the left
    const auto deck = deckBounds().toNearestInt();
    presetBar.setBounds (Rectangle<int> (360, 44).withCentre ({ kWidth / 2, deck.getY() + 40 }));
    abButton.setBounds (Rectangle<int> (deck.getRight() - 24 - 96, deck.getY() + 20, 96, 40).expanded (m));
    saveSlot.setBounds (Rectangle<int> (deck.getX() + 24, deck.getY() + 20, 96, 40).expanded (m));

    // Row 2: four main knobs, then four reserved smaller ones (75% size)
    const int bigKnob = 104, bigW = 120, bigH = bigKnob + 44;
    const int smallKnob = 78, smallW = 92, smallH = smallKnob + 44;
    const int rowTop = deck.getY() + 86;
    int cx = deck.getX() + 40 + bigW / 2;
    for (auto* k : { &volume, &dynamics, &hardness, &damper })
    {
        k->setBounds (cx - bigW / 2, rowTop, bigW, bigH);
        cx += bigW;
    }
    cx += 50 - bigW / 2 + smallW / 2;
    for (auto& s : smallKnobSlots)
    {
        s.setBounds (cx - smallW / 2, rowTop + (bigKnob - smallKnob) / 2, smallW, smallH);
        cx += smallW;
    }
    action.setBounds (smallKnobSlots[0].getBounds());
    age.setBounds (smallKnobSlots[1].getBounds());

    // Row 3: velocity slider, sustain and reserved buttons
    const int row3 = deck.getY() + 268;
    velocity.setBounds (deck.getX() + 30, row3, 400, 40);
    x = deck.getX() + 470;
    sustainButton.setBounds (Rectangle<int> (x, row3, 112, 40).expanded (m));
    x += 112 + 14;
    for (auto& s : buttonSlots)
    {
        s.setBounds (Rectangle<int> (x, row3, 112, 40).expanded (m));
        x += 112 + 14;
    }

    // Row 4: status
    status.setBounds (deck.getX() + 30, deck.getBottom() - 46, deck.getWidth() - 60, 30);
}
