# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024-2026 Tom Bleher, Igor Korover
#
# chargesharing_core: shared physics library for the LogA charge-sharing model,
# noise injection, and Gaussian position reconstruction. Consumed by both the
# standalone Geant4 simulation (epicChargeSharing) and the EICrecon plugin
# (chargeSharingRecon).

# Allow using this file as a sub-project (so the eicrecon plugin can
# add_subdirectory() it without having all of Geant4 in scope).
if(NOT DEFINED PROJECT_NAME)
    cmake_minimum_required(VERSION 3.16)
    project(chargesharing_core LANGUAGES CXX)
endif()

if(NOT ROOT_FOUND)
    find_package(ROOT CONFIG REQUIRED COMPONENTS Core MathCore Minuit Hist Matrix)
endif()

add_library(chargesharing_core STATIC
    src/ChargeSharingCore.cc
    src/NoiseModel.cc
    src/GaussianFit.cc)
add_library(chargesharing::core ALIAS chargesharing_core)

target_include_directories(chargesharing_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_compile_features(chargesharing_core PUBLIC cxx_std_20)

target_link_libraries(chargesharing_core PUBLIC
    ROOT::Core
    ROOT::MathCore
    ROOT::Minuit
    ROOT::Hist
    ROOT::Matrix)

set_target_properties(chargesharing_core PROPERTIES
    POSITION_INDEPENDENT_CODE ON)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(chargesharing_core PRIVATE -Wall -Wextra -Wpedantic)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(chargesharing_core PRIVATE -Wno-\#warnings)
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(chargesharing_core PRIVATE -Wno-cpp)
    endif()
endif()
