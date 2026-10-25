# The instrument engine: plain C++17, no JUCE, so it can be tested and benchmarked anywhere.
add_library(openpiano_engine STATIC
    include/openpiano/Instrument.h
    include/openpiano/ParallelFor.h
    include/openpiano/PianoInstrument.h
    src/ParallelFor.cpp
    src/piano/NoteTable.h
    src/piano/StringModel.h
    src/piano/StringModel.cpp
    src/piano/PianoInstrument.cpp)

target_include_directories(openpiano_engine PUBLIC include PRIVATE src)
set_target_properties(openpiano_engine PROPERTIES POSITION_INDEPENDENT_CODE ON)

find_package(Threads REQUIRED)
target_link_libraries(openpiano_engine PUBLIC Threads::Threads)

if(MSVC)
    target_compile_definitions(openpiano_engine PRIVATE _USE_MATH_DEFINES)
else()
    # No -ffast-math: the regression test relies on exact IEEE arithmetic.
    target_compile_options(openpiano_engine PRIVATE -Wall -Wextra)
endif()
