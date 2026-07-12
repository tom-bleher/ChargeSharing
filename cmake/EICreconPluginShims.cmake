# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2024-2026 Tom Bleher, Igor Korover
#
# Bridges this out-of-tree build to the plugin_* helper macros that ship with
# an installed EICrecon (<prefix>/lib/cmake/EICrecon/jana_plugin.cmake), so the
# per-detector CMakeLists can use the exact upstream plugin_add(),
# plugin_glob_all(), plugin_add_dd4hep(), plugin_add_acts(), and
# plugin_add_event_model() macros.
#
# find_package(EICrecon) must run before including this module; its package
# config already resolves JANA, podio, spdlog, and Microsoft.GSL, which the
# upstream macros rely on.

if(NOT EICrecon_FOUND)
    message(FATAL_ERROR
        "EICreconPluginShims requires find_package(EICrecon) to have succeeded first.")
endif()

find_file(EICRECON_JANA_PLUGIN_MODULE jana_plugin.cmake
    HINTS "${EICrecon_DIR}"
    NO_DEFAULT_PATH)

if(NOT EICRECON_JANA_PLUGIN_MODULE)
    message(FATAL_ERROR
        "jana_plugin.cmake not found next to EICreconConfig.cmake (${EICrecon_DIR}). "
        "Point CMAKE_PREFIX_PATH at a full EICrecon install whose lib/cmake/EICrecon "
        "directory contains jana_plugin.cmake (EICrecon >= 1.17; any recent build "
        "of https://github.com/eic/EICrecon installs it).")
endif()

# Upstream's top-level CMakeLists defines the install destinations consumed by
# plugin_add(). Out of tree we install into <prefix>/plugins so eicrecon finds
# the .so files via $EICrecon_MY/plugins or -Pjana:plugin_path.
if(NOT DEFINED PLUGIN_OUTPUT_DIRECTORY)
    set(PLUGIN_OUTPUT_DIRECTORY "plugins")
endif()
if(NOT DEFINED PLUGIN_LIBRARY_OUTPUT_DIRECTORY)
    set(PLUGIN_LIBRARY_OUTPUT_DIRECTORY "lib")
endif()

include("${EICRECON_JANA_PLUGIN_MODULE}")
