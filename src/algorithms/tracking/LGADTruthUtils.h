// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#pragma once

#include <edm4hep/MCParticle.h>

#include <cstddef>

namespace eicrecon {

/// Return the generated ancestor of a Geant4-created particle. A generated
/// particle is identified by a nonzero generator status; malformed ancestry is
/// bounded to avoid following a cyclic relation indefinitely.
inline edm4hep::MCParticle generatedAncestor(edm4hep::MCParticle particle) {
    constexpr std::size_t kMaxAncestryDepth = 1024;
    for (std::size_t depth = 0;
         particle.isAvailable() && particle.getGeneratorStatus() == 0 && particle.parents_size() > 0 &&
         depth < kMaxAncestryDepth;
         ++depth) {
        particle = particle.getParents(0);
    }
    return particle;
}

} // namespace eicrecon
