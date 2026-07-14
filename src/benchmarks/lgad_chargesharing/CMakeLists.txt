# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024-2026 Tom Bleher, Igor Korover

set(PLUGIN_NAME LGAD_chargesharing_benchmark)

plugin_add(${PLUGIN_NAME})
plugin_glob_all(${PLUGIN_NAME})
plugin_add_dd4hep(${PLUGIN_NAME})
plugin_add_event_model(${PLUGIN_NAME})
# Out-of-tree builds consume EICrecon service headers (services/log/...) from
# the installed include tree rather than the in-tree src/ directory.
plugin_include_directories(${PLUGIN_NAME} SYSTEM PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${EICrecon_INCLUDE_DIR}
    ${EICrecon_INCLUDE_DIR}/EICrecon
)
plugin_link_libraries(${PLUGIN_NAME}
    ROOT::Core
    ROOT::Tree
)

# Install the benchmark driver scripts next to the plugin so they can be
# invoked from $EICrecon_MY/bin/lgad_chargesharing.
install(PROGRAMS
    ${CMAKE_CURRENT_LIST_DIR}/run_b0.sh
    ${CMAKE_CURRENT_LIST_DIR}/validate.py
    DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/lgad_chargesharing
)
