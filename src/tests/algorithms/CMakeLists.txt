# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024-2026 Tom Bleher, Igor Korover
#
# Catch2-based unit tests for the LGAD charge-sharing algorithms. Runs inside
# eic-shell where Catch2 v3 is provided by the container. Tests link against
# the shared algorithm OBJECT library (lgad_chargesharing_algorithms) so they
# exercise the exact compiled sources used by the production plugins.

find_package(Catch2 3 QUIET)
if(NOT Catch2_FOUND)
    find_package(Catch2 QUIET)
endif()
if(NOT Catch2_FOUND)
    message(STATUS "Catch2 not found: skipping LGAD algorithm unit tests.")
    return()
endif()

include(Catch)

add_executable(lgad_chargesharing_tests
    test_LGADChargeSharingRecon.cc
    test_LGADGaussianClustering.cc
)

target_link_libraries(lgad_chargesharing_tests PRIVATE
    lgad_chargesharing::algorithms
    chargesharing::core
    Catch2::Catch2WithMain
)

target_include_directories(lgad_chargesharing_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/algorithms
)

catch_discover_tests(lgad_chargesharing_tests
    TEST_PREFIX "lgad_chargesharing::"
)
