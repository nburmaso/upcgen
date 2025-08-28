##########################################################################
# Copyright (C) 2021-2025, Nazar Burmasov, Evgeny Kryshen
#
# E-mail of the corresponding author: nazar.burmasov@cern.ch
#
# This file is a part of Upcgen
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
##########################################################################

cmake_minimum_required(VERSION 3.9)

set(CMAKE_CXX_STANDARD 17)

project(Upcgen LANGUAGES CXX C)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")

# debug build
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "-O0 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")

# release build
#set(CMAKE_CXX_FLAGS "-O2 -fpic")

option(BUILD_WITH_PYTHIA6 "Build with Pythia6 support" OFF)
option(BUILD_WITH_PYTHIA8 "Build with Pythia8 support" ON)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

find_package(ROOT REQUIRED)
if (ROOT_FOUND)
    message("Found ROOT " ${ROOT_VERSION} ":")
    message("ROOT includes: " ${ROOT_INCLUDE_DIRS})
    message("ROOT libraries: " ${ROOT_LIBRARIES})
    include_directories(${ROOT_INCLUDE_DIRS})
endif ()

find_package(GSL REQUIRED)
if (GSL_FOUND)
    message("Found GSL " ${GSL_VERSION} ":")
    message("GSL includes: " ${GSL_INCLUDE_DIRS})
    message("GSL libraries: " ${GSL_LIBRARIES})
    include_directories(${GSL_INCLUDE_DIRS})
endif ()

set(PYBIND11_FINDPYTHON ON)
add_subdirectory(extern/pybind11)

file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
include_directories(${PROJECT_SOURCE_DIR}/include)

add_definitions(-DCROSS_SEC_DIR=\"${CMAKE_SOURCE_DIR}/cross_sections\")

# note: root with pythia6 support is required for this option
if (BUILD_WITH_PYTHIA6)
    add_definitions(-DUSE_PYTHIA6)
    set(optLibs ${optLibs} $ENV{ROOTSYS}/lib/libEGPythia6.so)
endif ()

if (BUILD_WITH_PYTHIA8)
    add_definitions(-DUSE_PYTHIA8)
    find_package(Pythia8 REQUIRED)
    set(optLibs ${optLibs} ${PYTHIA8_LIBRARY})
    include_directories(${PYTHIA8_INCLUDE_DIR})
endif()

# static and ...
set(Upcgenlib "Upcgenlib")
add_library(${Upcgenlib} STATIC ${sources})

# ... shared library
set(UpcgenlibShared "UpcgenlibShared")
add_library(${UpcgenlibShared} SHARED ${sources})
target_link_libraries(${UpcgenlibShared} ${optLibs} ${ROOT_LIBRARIES} ${GSL_LIBRARIES})
set_target_properties(${UpcgenlibShared} PROPERTIES OUTPUT_NAME "Upcgenlib")

# executable
add_executable(upcgen main.cpp)
target_link_libraries(upcgen ${Upcgenlib} ${optLibs} ${ROOT_LIBRARIES} ${GSL_LIBRARIES} $ENV{ROOTSYS}/lib/libEG.so)

# executable
add_executable(test test.cpp)
target_link_libraries(test ${Upcgenlib} ${optLibs} ${ROOT_LIBRARIES} ${GSL_LIBRARIES} $ENV{ROOTSYS}/lib/libEG.so)
