# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024-2026 Tom Bleher, Igor Korover

get_filename_component(PLUGIN_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)
set(PLUGIN_NAME ${PLUGIN_NAME}_lgad_chargesharing)

plugin_add(${PLUGIN_NAME})
plugin_glob_all(${PLUGIN_NAME})
plugin_add_dd4hep(${PLUGIN_NAME})
plugin_add_acts(${PLUGIN_NAME})
plugin_add_event_model(${PLUGIN_NAME})
plugin_link_libraries(${PLUGIN_NAME}
    lgad_chargesharing::algorithms
    chargesharing::core
    # Upstream CKF/ActsToTracks/AmbiguitySolver implementations for the
    # parallel charge-sharing tracking chain.
    EICrecon::algorithms_tracking_library
)
