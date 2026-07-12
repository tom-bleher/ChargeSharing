// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#include <JANA/JApplication.h>

#include "LGADChargeSharingMonitor.h"

extern "C" {

void InitPlugin(JApplication* app) {
    InitJANAPlugin(app);
    app->Add(new eicrecon::LGADChargeSharingMonitor());
}

} // extern "C"
