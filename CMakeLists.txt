# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024-2026 Tom Bleher, Igor Korover
#
# Out-of-tree EICrecon plugin suite for AC-LGAD charge-sharing reconstruction
# on the B0 tracker. Models the upstream
# eic/EICrecon layout: a shared algorithm OBJECT library plus small
# per-detector plugin libraries that register their factories via InitPlugin().

cmake_minimum_required(VERSION 3.24)

project(epic_lgad_chargesharing LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_PREFIX_PATH "$ENV{JANA_HOME}/lib/JANA/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(EICrecon REQUIRED)
find_package(DD4hep   REQUIRED COMPONENTS DDCore DDRec)
find_package(ROOT     REQUIRED COMPONENTS Core MathCore Minuit Minuit2 Hist Tree Matrix)
find_package(Acts     REQUIRED COMPONENTS Core PluginDD4hep PluginJson)
find_package(EDM4HEP  REQUIRED)
find_package(EDM4EIC  REQUIRED)
find_package(podio    REQUIRED)
find_package(fmt      REQUIRED)

include(EICreconPluginShims)

# Shared physics library (vendored copy of epicChargeSharing's core/).
if(NOT TARGET chargesharing::core)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/core
                     ${CMAKE_CURRENT_BINARY_DIR}/core)
endif()

add_subdirectory(src/algorithms)
add_subdirectory(src/detectors/B0TRK)
add_subdirectory(src/benchmarks/lgad_chargesharing)

option(BUILD_TESTING "Build Catch2 unit tests" ON)
if(BUILD_TESTING)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/tests/algorithms/CMakeLists.txt")
        enable_testing()
        add_subdirectory(src/tests/algorithms)
    else()
        message(STATUS "BUILD_TESTING=ON but src/tests/algorithms/ not present; skipping.")
    endif()
endif()
