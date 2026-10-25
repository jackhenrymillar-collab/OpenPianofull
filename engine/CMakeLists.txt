# The instrument engine: plain C++17, no JUCE, so it can be tested and benchmarked anywhere.
add_library(openpiano_engine STATIC
    include/openpiano/Instrument.h
    include/openpiano/ParallelFor.h
    include/openpiano/LegacyPiano.h
    include/openpiano/PianoControls.h
    include/openpiano/PianoInstrument.h
    src/ParallelFor.cpp
    src/piano/NoteTable.h
    src/piano/StringModel.h
    src/piano/StringModel.cpp
    src/piano/LegacyPiano.cpp
    src/piano/ScaleDesign.h
    src/piano/ScaleDesign.cpp
    src/piano/ModalPiano.h
    src/piano/ModalPiano.cpp
    src/piano/PianoInstrument.cpp)

target_include_directories(openpiano_engine PUBLIC include PRIVATE src)
set_target_properties(openpiano_engine PROPERTIES POSITION_INDEPENDENT_CODE ON)

find_package(Threads REQUIRED)
target_link_libraries(openpiano_engine PUBLIC Threads::Threads)

if(MSVC)
    target_compile_definitions(openpiano_engine PRIVATE _USE_MATH_DEFINES)
else()
    # No -ffast-math on the legacy model: its regression test relies on exact IEEE arithmetic.
    target_compile_options(openpiano_engine PRIVATE -Wall -Wextra)
    # The modal piano's inner loops vectorise much better with relaxed float rules.
    set_source_files_properties(src/piano/ModalPiano.cpp PROPERTIES COMPILE_OPTIONS "-O3;-ffast-math;-fno-finite-math-only")
endif()
