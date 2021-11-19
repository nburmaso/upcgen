cmake_minimum_required(VERSION 2.8)
project(LepGenerator)

set(CMAKE_CXX_STANDARD 11)

option(BUILD_WITH_PYTHIA6 "Build with Pythia6 support" OFF)
option(BUILD_WITH_PYTHIA8 "Build with Pythia8 support" OFF)
option(BUILD_WITH_OPENMP "Build with OpenMP support" OFF)

find_package(OpenMP)
if (OPENMP_FOUND AND BUILD_WITH_OPENMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS}
        $ENV{PYTHIA8}/include)

include_directories(${PROJECT_SOURCE_DIR}/include
        ${ROOT_INCLUDE_DIRS})

file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_executable(generate main.cpp ${sources} ${headers})

target_link_libraries(generate
        ${ROOT_LIBRARIES}
        $ENV{ROOTSYS}/lib/libEG.so)

if (BUILD_WITH_PYTHIA6)
    add_definitions(-DUSE_PYTHIA6)
    target_link_libraries(generate $ENV{ROOTSYS}/lib/libEGPythia6.so)
endif ()

if (BUILD_WITH_PYTHIA8)
    add_definitions(-DUSE_PYTHIA8)
    target_link_libraries(generate
            $ENV{PYTHIA8}/lib/libpythia8.so
            $ENV{ROOTSYS}/lib/libEGPythia8.so)
endif()