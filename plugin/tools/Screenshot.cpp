/*
OpenPiano: an open source piano engine based on physical modeling
Copyright (C) 2021-2022 Michele Perrone, and contributors
Licensed under the GNU Affero General Public License v3 or later (see LICENSE).
*/

// Saves a PNG of the plugin window without a host: openpiano-screenshot out.png [scale]
// Useful for reviewing UI changes on a build machine. On Linux, run it under xvfb-run.

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_gui_basics/juce_gui_basics.h>


int main (int argc, char** argv)
{
    juce::ScopedJuceInitialiser_GUI gui;
    const juce::File out (juce::File::getCurrentWorkingDirectory().getChildFile (argc > 1 ? argv[1] : "screenshot.png"));
    const float scale = argc > 2 ? (float) std::atof (argv[2]) : 2.0f;

    std::unique_ptr<juce::AudioProcessor> processor (createPluginFilter());
    processor->setRateAndBufferSizeDetails (48000.0, 256);
    processor->prepareToPlay (48000.0, 256);

    std::unique_ptr<juce::AudioProcessorEditor> editor (processor->createEditor());
    editor->setVisible (true);

    juce::Image image (juce::Image::ARGB, juce::roundToInt (editor->getWidth() * scale),
                       juce::roundToInt (editor->getHeight() * scale), true);
    {
        juce::Graphics g (image);
        g.addTransform (juce::AffineTransform::scale (scale));
        editor->paintEntireComponent (g, true);
    }

    out.deleteFile();
    juce::FileOutputStream stream (out);
    juce::PNGImageFormat png;
    if (! stream.openedOk() || ! png.writeImageToStream (image, stream))
        return 1;

    std::printf ("Wrote %s\n", out.getFullPathName().toRawUTF8());
    editor = nullptr;
    processor->releaseResources();
    return 0;
}
