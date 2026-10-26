cmake_minimum_required(VERSION 3.22)

project(OpenPianoFull VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(APPLE)
    # Universal binary (Apple silicon + Intel), macOS 11 and newer
    set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0" CACHE STRING "Minimum macOS version")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(OPENPIANO_BUILD_PLUGIN "Build the JUCE plugin (AU, VST3, Standalone)" ON)
option(OPENPIANO_BUILD_TOOLS "Build the render tool and regression test" ON)

add_subdirectory(engine)

if(OPENPIANO_BUILD_TOOLS)
    enable_testing()
    add_subdirectory(tools)
endif()

if(OPENPIANO_BUILD_PLUGIN)
    add_subdirectory(plugin)
endif()
