# The JUCE plugin. JUCE is downloaded at configure time. To use a local copy instead, pass
#   -DFETCHCONTENT_SOURCE_DIR_JUCE=/path/to/JUCE
include(FetchContent)
FetchContent_Declare(JUCE
    GIT_REPOSITORY https://github.com/juce-framework/JUCE.git
    GIT_TAG 8.0.9
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(JUCE)

# Plugin identity. These are placeholders until the final product name is chosen.
# Changing PLUGIN_CODE or MANUFACTURER_CODE makes hosts treat it as a different plugin.
set(OPENPIANO_PRODUCT_NAME "Nocturne" CACHE STRING "Plugin name shown in hosts")
set(OPENPIANO_COMPANY_NAME "Openpianofull" CACHE STRING "Company name shown in hosts")
set(OPENPIANO_MANUFACTURER_CODE "Opfl" CACHE STRING "Four-character manufacturer code")
set(OPENPIANO_PLUGIN_CODE "Opfp" CACHE STRING "Four-character plugin code, exactly one upper-case letter")

set(OPENPIANO_FORMATS VST3 Standalone)
if(APPLE)
    list(APPEND OPENPIANO_FORMATS AU)
endif()

juce_add_plugin(OpenPianoPlugin
    PRODUCT_NAME "${OPENPIANO_PRODUCT_NAME}"
    COMPANY_NAME "${OPENPIANO_COMPANY_NAME}"
    BUNDLE_ID "com.openpianofull.nocturne"
    PLUGIN_MANUFACTURER_CODE ${OPENPIANO_MANUFACTURER_CODE}
    PLUGIN_CODE ${OPENPIANO_PLUGIN_CODE}
    IS_SYNTH TRUE
    NEEDS_MIDI_INPUT TRUE
    NEEDS_MIDI_OUTPUT FALSE
    IS_MIDI_EFFECT FALSE
    EDITOR_WANTS_KEYBOARD_FOCUS FALSE
    COPY_PLUGIN_AFTER_BUILD FALSE
    AU_MAIN_TYPE kAudioUnitType_MusicDevice
    VST3_CATEGORIES Instrument Synth
    MICROPHONE_PERMISSION_ENABLED FALSE
    FORMATS ${OPENPIANO_FORMATS})

juce_add_binary_data(OpenPianoAssets
    NAMESPACE OpenPianoAssets
    SOURCES
        Resources/fonts/CormorantGaramond.ttf
        Resources/fonts/Jost.ttf)

target_sources(OpenPianoPlugin PRIVATE
    Source/Parameters.h
    Source/Parameters.cpp
    Source/PluginProcessor.h
    Source/PluginProcessor.cpp
    Source/PluginEditor.h
    Source/PluginEditor.cpp
    Source/ui/Theme.h
    Source/ui/Theme.cpp
    Source/ui/Sky.h
    Source/ui/Sky.cpp
    Source/ui/Widgets.h
    Source/ui/Widgets.cpp)

target_include_directories(OpenPianoPlugin PRIVATE Source)

target_compile_definitions(OpenPianoPlugin PUBLIC
    JUCE_WEB_BROWSER=0
    JUCE_USE_CURL=0
    JUCE_VST3_CAN_REPLACE_VST2=0
    JUCE_DISPLAY_SPLASH_SCREEN=0)

target_link_libraries(OpenPianoPlugin
    PRIVATE
        openpiano_engine
        OpenPianoAssets
        juce::juce_audio_utils
        juce::juce_audio_processors
    PUBLIC
        juce::juce_recommended_config_flags
        juce::juce_recommended_lto_flags
        juce::juce_recommended_warning_flags)

# Optional: a command-line tool that saves a PNG of the plugin window
option(OPENPIANO_BUILD_SCREENSHOT "Build the window screenshot tool" OFF)
if(OPENPIANO_BUILD_SCREENSHOT)
    juce_add_console_app(OpenPianoScreenshot PRODUCT_NAME "openpiano-screenshot")
    target_sources(OpenPianoScreenshot PRIVATE tools/Screenshot.cpp)
    target_include_directories(OpenPianoScreenshot PRIVATE Source)
    target_compile_definitions(OpenPianoScreenshot PRIVATE JUCE_WEB_BROWSER=0 JUCE_USE_CURL=0)
    target_link_libraries(OpenPianoScreenshot PRIVATE OpenPianoPlugin openpiano_engine juce::juce_gui_basics juce::juce_audio_utils juce::juce_recommended_config_flags)
endif()
