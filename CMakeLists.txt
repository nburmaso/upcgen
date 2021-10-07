cmake_minimum_required(VERSION 2.8)
project(LepGenerator)

set(CMAKE_CXX_STANDARD 11)

option(BUILD_WITH_PYTHIA6 "Build with Pythia6 support" ON)

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS}
        $ENV{PYTHIA8}/include)

# use in case if GSL is not fetched from system
# set (GSL_ROOT_DIR "path/to/install")
find_package(GSL REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include
        ${ROOT_INCLUDE_DIRS})

file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_executable(generate main.cpp ${sources} ${headers})

target_link_libraries(generate
        ${ROOT_LIBRARIES}
        $ENV{ROOTSYS}/lib/libEG.so
        $ENV{PYTHIA8}/lib/libpythia8.so
        $ENV{ROOTSYS}/lib/libEGPythia8.so
        GSL::gsl)

if (BUILD_WITH_PYTHIA6)
    add_definitions(-DUSE_PYTHIA6)
    target_link_libraries(generate $ENV{ROOTSYS}/lib/libEGPythia6.so)
endif ()