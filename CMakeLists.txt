cmake_minimum_required(VERSION 2.8)
project(LepGenerator)

set(CMAKE_CXX_STANDARD 11)

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS}
        $ENV{PYTHIA8}/include)

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
        GSL::gsl)