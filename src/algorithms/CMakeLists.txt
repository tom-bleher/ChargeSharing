# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024-2026 Tom Bleher, Igor Korover
#
# Shared OBJECT library holding the LGAD charge-sharing algorithms.
# Linked into each per-detector plugin library (under src/detectors/*) and
# into the unit test executable.

add_library(lgad_chargesharing_algorithms OBJECT
    reco/LGADChargeSharingRecon.cc
    tracking/LGADGaussianClustering.cc
)

target_include_directories(lgad_chargesharing_algorithms PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${EICrecon_INCLUDE_DIR}
    ${EICrecon_INCLUDE_DIR}/EICrecon
)

target_link_libraries(lgad_chargesharing_algorithms PUBLIC
    chargesharing::core
    algorithms::core
    EDM4HEP::edm4hep
    EDM4EIC::edm4eic
    podio::podio
    DD4hep::DDCore
    DD4hep::DDRec
    Acts::Core
    ROOT::Core
    ROOT::MathCore
    ROOT::Minuit2
    ROOT::Hist
    ROOT::Matrix
    fmt::fmt
)

target_compile_definitions(lgad_chargesharing_algorithms PUBLIC
    JANA_VERSION_MAJOR=${JANA_VERSION_MAJOR}
    JANA_VERSION_MINOR=${JANA_VERSION_MINOR}
    JANA_VERSION_PATCH=${JANA_VERSION_PATCH}
    EICRECON_PLUGIN_BUILD=1
)

set_target_properties(lgad_chargesharing_algorithms PROPERTIES
    POSITION_INDEPENDENT_CODE ON
)

add_library(lgad_chargesharing::algorithms ALIAS lgad_chargesharing_algorithms)
