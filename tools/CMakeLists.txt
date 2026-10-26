# dr_wav, for writing renders
add_library(openpiano_drwav STATIC ../third_party/dr_wav/dr_wav.cpp)
target_include_directories(openpiano_drwav PUBLIC ../third_party/dr_wav)

# The original OpenPiano code, frozen as a reference for regression tests
add_library(openpiano_upstream STATIC upstream/array_helpers.cpp)
target_include_directories(openpiano_upstream PUBLIC upstream)
target_link_libraries(openpiano_upstream PUBLIC openpiano_drwav)

add_library(openpiano_harness STATIC
    harness/Harness.h
    harness/Harness.cpp
    harness/UpstreamRender.cpp)
target_include_directories(openpiano_harness PUBLIC harness)
target_link_libraries(openpiano_harness PUBLIC openpiano_engine openpiano_upstream openpiano_drwav)

add_executable(openpiano-render harness/RenderMain.cpp)
target_link_libraries(openpiano-render PRIVATE openpiano_harness)
# For --scale, which prints the piano's design table
target_include_directories(openpiano-render PRIVATE ../engine/src)

add_executable(openpiano-regression-test ../tests/RegressionTest.cpp)
target_link_libraries(openpiano-regression-test PRIVATE openpiano_harness)
add_test(NAME regression_matches_upstream COMMAND openpiano-regression-test)

add_executable(openpiano-piano-test ../tests/PianoTest.cpp)
target_link_libraries(openpiano-piano-test PRIVATE openpiano_engine)
target_include_directories(openpiano-piano-test PRIVATE ../engine/src)
add_test(NAME piano_checks COMMAND openpiano-piano-test)

# Single held notes, for calibrating against recorded samples (tools/calibrate)
add_executable(openpiano-notes harness/NotesMain.cpp)
target_link_libraries(openpiano-notes PRIVATE openpiano_harness)

# Renders a text list of note and pedal events (for example a converted MIDI file)
add_executable(openpiano-events harness/EventsMain.cpp)
target_link_libraries(openpiano-events PRIVATE openpiano_harness)
