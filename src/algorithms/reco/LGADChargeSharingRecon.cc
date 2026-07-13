// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#include "LGADChargeSharingRecon.h"

#include <DD4hep/Detector.h>
#include <DD4hep/Alignments.h>
#include <DD4hep/Readout.h>
#include <DD4hep/Segmentations.h>
#include <DD4hep/Objects.h>
#include <DD4hep/Volumes.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DDSegmentation/BitFieldCoder.h>
#include <DDSegmentation/CartesianGridXY.h>
#include <DDSegmentation/CartesianGridXZ.h>
#include <DDSegmentation/Segmentation.h>
#include <Evaluator/DD4hepUnits.h>
#include <algorithms/geo.h>
#include <edm4eic/CovDiag3f.h>
#include <edm4hep/Vector3f.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr double kElementaryChargeC = 1.602176634e-19;

/// Detector-specific DD4hep constants. Used only when segmentation does not
/// supply the value directly (e.g. sensor thickness). Must stay in sync with
/// the compact XML under `eic/epic`.
struct DetectorConstants {
    std::string siliconThickness;
    std::string detectorSize;
    std::string pixelSize;
    std::string copperThickness;
};

DetectorConstants getDetectorConstants(const std::string& readout) {
    DetectorConstants constants;
    if (readout.find("LumiSpec") != std::string::npos || readout.find("Lumi") != std::string::npos) {
        constants.siliconThickness = "LumiSpecTracker_Si_DZ";
        constants.detectorSize = "LumiSpecTracker_DXY";
        constants.pixelSize = "LumiSpecTracker_pixelSize";
        constants.copperThickness = "LumiSpecTracker_Cu_DZ";
    } else if (readout.find("B0") != std::string::npos) {
        constants.siliconThickness = "B0TrackerSensorThickness";
        constants.pixelSize = "B0TrackerElectrodeSize";
    }
    return constants;
}

double getSensorThicknessFromReadout(const dd4hep::Detector* detector, const std::string& readoutName) {
    if (!detector)
        return 0.0;
    try {
        dd4hep::Readout readout = detector->readout(readoutName);
        if (!readout.isValid())
            return 0.0;
        dd4hep::IDDescriptor idDesc = readout.idSpec();
        if (!idDesc.isValid())
            return 0.0;

        std::string detName = readoutName;
        if (detName.size() > 4 && detName.substr(detName.size() - 4) == "Hits") {
            detName = detName.substr(0, detName.size() - 4);
        }

        dd4hep::DetElement det = detector->detector(detName);
        if (!det.isValid()) {
            det = detector->detector(readoutName);
        }

        if (det.isValid() && det.placement().isValid()) {
            dd4hep::Volume vol = det.placement().volume();
            if (vol.isValid() && vol.solid().isValid()) {
                auto dims = vol.solid().dimensions();
                if (dims.size() >= 3) {
                    return 2.0 * dims[2] / dd4hep::mm;
                }
            }
        }
    } catch (...) {
    }
    return 0.0;
}

void applyDetectorFallbacks(const std::string& readout, eicrecon::LGADChargeSharingRecon::Geometry& geom) {
    if (readout.find("B0") != std::string::npos) {
        geom.detectorThicknessMM = 0.3;
    }
}

} // namespace

namespace eicrecon {

namespace core = ::chargesharing::core;

// ---------------------------------------------------------------------------
// init()
// ---------------------------------------------------------------------------
void LGADChargeSharingRecon::init() {
    if (m_cfg.neighborhoodRadius < 0) {
        m_cfg.neighborhoodRadius = 0;
    }

    // Noise model setup
    core::NoiseConfig noiseConfig;
    noiseConfig.enabled = m_cfg.noiseEnabled;
    noiseConfig.gainSigmaMin = 0.01;
    noiseConfig.gainSigmaMax = 0.05;
    noiseConfig.electronNoiseCount = m_cfg.noiseElectronCount;
    noiseConfig.elementaryCharge = kElementaryChargeC;
    m_noise_model.setConfig(noiseConfig);

    if (m_skip_dd4hep_init) {
        info("LGADChargeSharingRecon: DD4hep init skipped (test geometry injected)");
        return;
    }

    // ------------------------------------------------------------------
    // DD4hep geometry lookup via algorithms::GeoSvc singleton.
    // Pitch, offset and index fields come from readout segmentation. Physical
    // electrode size and sensor thickness come from compact-XML constants.
    // ------------------------------------------------------------------
    const auto& geo = algorithms::GeoSvc::instance();
    const dd4hep::Detector* detector = geo.detector();
    m_converter = geo.cellIDPositionConverter();

    if (!detector) {
        warning("DD4hep detector unavailable; LGADChargeSharingRecon will run with defaults");
        return;
    }

    try {
        dd4hep::Readout readout = detector->readout(m_cfg.readout);
        dd4hep::Segmentation segmentation = readout.segmentation();

        if (!segmentation.isValid()) {
            warning("Readout '{}' has no valid segmentation", m_cfg.readout);
        } else {
            const auto* segImplXY =
                dynamic_cast<const dd4hep::DDSegmentation::CartesianGridXY*>(segmentation.segmentation());
            const dd4hep::DDSegmentation::CartesianGridXZ* segImplXZ = nullptr;
            if (segImplXY == nullptr) {
                segImplXZ =
                    dynamic_cast<const dd4hep::DDSegmentation::CartesianGridXZ*>(segmentation.segmentation());
                m_geom.useXZCoordinates = (segImplXZ != nullptr);
            }

            if (segImplXY == nullptr && segImplXZ == nullptr) {
                warning("Segmentation for readout '{}' is '{}'; expected CartesianGridXY or CartesianGridXZ",
                        m_cfg.readout, segmentation.type());
            } else {
                m_decoder = segmentation.decoder();
                m_segmentation = segmentation.segmentation();

                if (segImplXY != nullptr) {
                    m_geom.pixelSpacingXMM = segImplXY->gridSizeX() / dd4hep::mm;
                    m_geom.pixelSpacingYMM = segImplXY->gridSizeY() / dd4hep::mm;
                    m_geom.gridOffsetXMM = segImplXY->offsetX() / dd4hep::mm;
                    m_geom.gridOffsetYMM = segImplXY->offsetY() / dd4hep::mm;
                    m_geom.fieldNameX = segImplXY->fieldNameX();
                    m_geom.fieldNameY = segImplXY->fieldNameY();
                } else {
                    m_geom.pixelSpacingXMM = segImplXZ->gridSizeX() / dd4hep::mm;
                    m_geom.pixelSpacingYMM = segImplXZ->gridSizeZ() / dd4hep::mm;
                    m_geom.gridOffsetXMM = segImplXZ->offsetX() / dd4hep::mm;
                    m_geom.gridOffsetYMM = segImplXZ->offsetZ() / dd4hep::mm;
                    m_geom.fieldNameX = segImplXZ->fieldNameX();
                    m_geom.fieldNameY = segImplXZ->fieldNameZ();
                }

                info("Using DD4hep {} segmentation for readout '{}': pitch=({}, {}) mm, "
                     "offset=({}, {}) mm",
                     m_geom.useXZCoordinates ? "CartesianGridXZ" : "CartesianGridXY", m_cfg.readout,
                     m_geom.pixelSpacingXMM, m_geom.pixelSpacingYMM, m_geom.gridOffsetXMM,
                     m_geom.gridOffsetYMM);
            }
        }

        // Optional per-detector constants for sensor thickness / pad size overrides.
        const auto detConstants = getDetectorConstants(m_cfg.readout);
        if (!detConstants.siliconThickness.empty()) {
            try {
                m_geom.detectorThicknessMM =
                    detector->constantAsDouble(detConstants.siliconThickness) / dd4hep::mm;
            } catch (const std::exception& ex) {
                debug("Optional constant '{}' not available ({}); keeping default thickness.",
                      detConstants.siliconThickness, ex.what());
            }
        }
        if (!detConstants.pixelSize.empty()) {
            try {
                double p = detector->constantAsDouble(detConstants.pixelSize) / dd4hep::mm;
                if (p > 0.0) {
                    m_geom.pixelSizeXMM = p;
                    m_geom.pixelSizeYMM = p;
                }
            } catch (const std::exception& ex) {
                warning("Electrode-size constant '{}' not available ({}); keeping model fallback {} mm.",
                        detConstants.pixelSize, ex.what(), m_geom.pixelSizeXMM);
            }
        }
        if (!detConstants.copperThickness.empty()) {
            try {
                m_geom.pixelThicknessMM =
                    detector->constantAsDouble(detConstants.copperThickness) / dd4hep::mm;
            } catch (const std::exception& ex) {
                debug("Optional constant '{}' not available ({}); keeping default pixel thickness.",
                      detConstants.copperThickness, ex.what());
            }
        }

        if (detConstants.siliconThickness.empty() && detConstants.pixelSize.empty()) {
            double volumeThickness = getSensorThicknessFromReadout(detector, m_cfg.readout);
            if (volumeThickness > 0.0 && volumeThickness < 10.0) {
                m_geom.detectorThicknessMM = volumeThickness;
            } else {
                applyDetectorFallbacks(m_cfg.readout, m_geom);
                debug("Sensor thickness not derivable from DD4hep volume (got {} mm); "
                      "using hard-coded fallback {} mm for readout '{}'.",
                      volumeThickness, m_geom.detectorThicknessMM, m_cfg.readout);
            }
        }
    } catch (const std::exception& ex) {
        warning("Failed to derive segmentation for readout '{}': {}", m_cfg.readout, ex.what());
    }

    info("Final geometry: pitch=({}, {}) mm, pad=({}, {}) mm, Si thickness={} mm",
         m_geom.pixelSpacingXMM, m_geom.pixelSpacingYMM, m_geom.pixelSizeXMM, m_geom.pixelSizeYMM,
         m_geom.detectorThicknessMM);
}

std::pair<int, int> LGADChargeSharingRecon::indexRangeForBox(double halfExtentMM,
                                                              double pitchMM,
                                                              double offsetMM) {
    if (!(std::isfinite(halfExtentMM) && halfExtentMM > 0.0 && std::isfinite(pitchMM) &&
          pitchMM > 0.0 && std::isfinite(offsetMM))) {
        throw std::invalid_argument("Invalid sensor extent, pitch, or offset");
    }

    const double lower = std::ceil((-halfExtentMM - offsetMM) / pitchMM);
    const double upper = std::floor((halfExtentMM - offsetMM) / pitchMM);
    if (lower < static_cast<double>(std::numeric_limits<int>::min()) ||
        upper > static_cast<double>(std::numeric_limits<int>::max()) || lower > upper) {
        throw std::out_of_range("Sensor extent does not define a valid integer cell range");
    }
    return {static_cast<int>(lower), static_cast<int>(upper)};
}

// ---------------------------------------------------------------------------
// process()
// ---------------------------------------------------------------------------
void LGADChargeSharingRecon::process(const Input& input, const Output& output) const {
    const auto [sim_hits] = input;
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
    auto [raw_hits, rec_hits, links, assocs] = output;
#else
    auto [raw_hits, rec_hits, assocs] = output;
#endif

    const bool useXZ = m_geom.useXZCoordinates;
    const double pitchX = m_geom.pixelSpacingXMM;
    const double pitchY = m_geom.pixelSpacingYMM;
    const double defaultVarX = (pitchX * pitchX) / 12.0;
    const double defaultVarY = (pitchY * pitchY) / 12.0;
    const double sensorT = m_geom.detectorThicknessMM;
    const double varZ = (sensorT * sensorT) / 12.0;

    edm4eic::CovDiag3f posError;
    if (useXZ) {
        posError = edm4eic::CovDiag3f{static_cast<float>(defaultVarX), static_cast<float>(varZ),
                                      static_cast<float>(defaultVarY)};
    } else {
        posError = edm4eic::CovDiag3f{static_cast<float>(defaultVarX), static_cast<float>(defaultVarY),
                                      static_cast<float>(varZ)};
    }

    // ------------------------------------------------------------------
    // Pass 1: spread every SimTrackerHit over its pad neighborhood and record
    // the noiseless induced charge on each pad. No electronics yet.
    // ------------------------------------------------------------------
    std::vector<PadContribution> contributions;
    contributions.reserve(sim_hits->size());

    for (std::size_t hitIndex = 0; hitIndex < sim_hits->size(); ++hitIndex) {
        const auto hit = (*sim_hits)[hitIndex];
        const double edep = hit.getEDep();
        if (edep < static_cast<double>(m_cfg.minEDepGeV)) {
            continue;
        }

        SingleHitInput singleInput{};
        const auto pos = hit.getPosition();
        if (m_converter == nullptr) {
            warning("CellID position converter unavailable; dropping hit {:#018x}", hit.getCellID());
            continue;
        }
        const auto* context = m_converter->findContext(hit.getCellID());
        if (context == nullptr) {
            warning("No DetElement context for hit {:#018x}; dropping it", hit.getCellID());
            continue;
        }
        const auto localPos = context->element.nominal().worldToLocal(
            dd4hep::Position(pos.x * dd4hep::mm, pos.y * dd4hep::mm, pos.z * dd4hep::mm));

        try {
            dd4hep::Box sensorBox = context->element.solid();
            const auto xRange = indexRangeForBox(sensorBox->GetDX() / dd4hep::mm,
                                                  m_geom.pixelSpacingXMM,
                                                  m_geom.gridOffsetXMM);
            const double secondHalfExtent =
                (useXZ ? sensorBox->GetDZ() : sensorBox->GetDY()) / dd4hep::mm;
            const auto yRange = indexRangeForBox(secondHalfExtent, m_geom.pixelSpacingYMM,
                                                  m_geom.gridOffsetYMM);
            singleInput.indexBounds = SingleHitInput::IndexBounds{
                xRange.first, xRange.second, yRange.first, yRange.second};
        } catch (const std::exception& ex) {
            error("Sensitive DetElement for hit {:#018x} is not a supported Box: {}; dropping hit",
                  hit.getCellID(), ex.what());
            continue;
        }

        std::uint64_t centerCellID = hit.getCellID();
        if (m_segmentation != nullptr) {
            centerCellID = m_segmentation->cellID(localPos, dd4hep::Position{}, centerCellID);
        }
        if (useXZ) {
            // CartesianGridXZ: map (x, y, z) -> local (x, z, y) so that the
            // algorithm's local Y is the detector's Z coordinate.
            singleInput.hitPositionMM = {localPos.x() / dd4hep::mm, localPos.z() / dd4hep::mm,
                                         localPos.y() / dd4hep::mm};
        } else {
            singleInput.hitPositionMM = {localPos.x() / dd4hep::mm, localPos.y() / dd4hep::mm,
                                         localPos.z() / dd4hep::mm};
        }
        singleInput.energyDepositGeV = edep;
        singleInput.cellID = centerCellID;

        if (m_decoder != nullptr) {
            const int idxI = static_cast<int>(m_decoder->get(centerCellID, m_geom.fieldNameX));
            const int idxJ = static_cast<int>(m_decoder->get(centerCellID, m_geom.fieldNameY));
            singleInput.pixelIndexHint = std::pair<int, int>{idxI, idxJ};
            const auto globalCenter = m_converter->position(centerCellID);
            const auto center = context->element.nominal().worldToLocal(globalCenter);
            if (useXZ) {
                singleInput.pixelHintMM = std::array<double, 3>{center.x() / dd4hep::mm,
                                                                 center.z() / dd4hep::mm,
                                                                 center.y() / dd4hep::mm};
            } else {
                singleInput.pixelHintMM = std::array<double, 3>{center.x() / dd4hep::mm,
                                                                 center.y() / dd4hep::mm,
                                                                 center.z() / dd4hep::mm};
            }
        } else if (m_converter != nullptr) {
            const auto globalCenter = m_converter->position(hit.getCellID());
            const auto center = context->element.nominal().worldToLocal(globalCenter);
            if (useXZ) {
                singleInput.pixelHintMM = std::array<double, 3>{center.x() / dd4hep::mm,
                                                                 center.z() / dd4hep::mm,
                                                                 center.y() / dd4hep::mm};
            } else {
                singleInput.pixelHintMM = std::array<double, 3>{center.x() / dd4hep::mm,
                                                                 center.y() / dd4hep::mm,
                                                                 center.z() / dd4hep::mm};
            }
        }

        const auto result = processSingleHitImpl(singleInput);
        if (!(result.totalCollectedChargeC > 0.0)) {
            continue;
        }

        for (const auto& neighbor : result.neighbors) {
            if (!(neighbor.chargeC > 0.0)) {
                continue;
            }

            std::uint64_t padCellID = centerCellID;
            if (m_decoder != nullptr) {
                m_decoder->set(padCellID, m_geom.fieldNameX, result.pixelRowIndex + neighbor.di);
                m_decoder->set(padCellID, m_geom.fieldNameY, result.pixelColIndex + neighbor.dj);
            } else if (neighbor.di != 0 || neighbor.dj != 0) {
                // A decoder is required to assign distinct neighbouring cells.
                continue;
            }

            edm4hep::Vector3f padPosition;
            if (m_converter != nullptr) {
                const auto globalCenter = m_converter->position(padCellID);
                padPosition = edm4hep::Vector3f{static_cast<float>(globalCenter.x() / dd4hep::mm),
                                                static_cast<float>(globalCenter.y() / dd4hep::mm),
                                                static_cast<float>(globalCenter.z() / dd4hep::mm)};
            } else if (useXZ) {
                padPosition = edm4hep::Vector3f{static_cast<float>(neighbor.pixelXMM),
                                                static_cast<float>(result.nearestPixelCenterMM[2]),
                                                static_cast<float>(neighbor.pixelYMM)};
            } else {
                padPosition = edm4hep::Vector3f{static_cast<float>(neighbor.pixelXMM),
                                                static_cast<float>(neighbor.pixelYMM),
                                                static_cast<float>(result.nearestPixelCenterMM[2])};
            }

            PadContribution contribution{};
            contribution.cellID = padCellID;
            contribution.position = padPosition;
            contribution.timeNs = hit.getTime();
            contribution.inducedChargeC = neighbor.chargeC;
            contribution.inducedEnergyGeV = edep * neighbor.fraction;
            contribution.simHitIndex = static_cast<int>(hitIndex);
            contributions.push_back(contribution);
        }
    }

    // ------------------------------------------------------------------
    // Pass 2: aggregate contributions into readout channels (sum before
    // electronics), then apply noise, threshold and digitization once per
    // channel and emit the raw/reconstructed hits.
    // ------------------------------------------------------------------
    const auto channels = aggregateChannels(contributions, m_cfg.integrationWindowNs);
    const double thresholdChargeC = m_cfg.thresholdElectrons * kElementaryChargeC;

    for (const auto& channel : channels) {
        double chargeC = channel.chargeC;
        if (m_cfg.noiseEnabled) {
            chargeC = m_noise_model.applyNoise(chargeC);
        }
        if (!(chargeC > thresholdChargeC)) {
            continue;
        }

        const double chargeElectrons = chargeC / kElementaryChargeC;
        const auto rawCharge = static_cast<std::int32_t>(std::clamp(
            std::llround(chargeElectrons),
            static_cast<long long>(std::numeric_limits<std::int32_t>::min()),
            static_cast<long long>(std::numeric_limits<std::int32_t>::max())));

        auto raw = raw_hits->create(channel.cellID, rawCharge,
                                    static_cast<std::int32_t>(channel.timeNs * 1e3));
        auto rec = rec_hits->create(channel.cellID, channel.position, posError,
                                    static_cast<float>(channel.timeNs), 0.0f,
                                    static_cast<float>(channel.energyGeV), 0.0f);
        rec.setRawHit(raw);

        for (const auto& [simHitIndex, weight] : channel.contributors) {
            const auto simHit = (*sim_hits)[static_cast<std::size_t>(simHitIndex)];
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
            auto link = links->create();
            link.setFrom(raw);
            link.setTo(simHit);
            link.setWeight(static_cast<float>(weight));
#endif
            auto assoc = assocs->create();
            assoc.setRawHit(raw);
            assoc.setSimHit(simHit);
            assoc.setWeight(static_cast<float>(weight));
        }
    }
}

// ---------------------------------------------------------------------------
// Channel aggregation (pure; unit-tested).
// ---------------------------------------------------------------------------
std::vector<LGADChargeSharingRecon::AggregatedChannel>
LGADChargeSharingRecon::aggregateChannels(const std::vector<PadContribution>& contributions,
                                          double integrationWindowNs) {
    // Accumulator carrying the extra running sums that the public channel type
    // does not need to expose (time anchor, charge-weighted time numerator).
    struct Bucket {
        AggregatedChannel channel;
        double anchorTimeNs{0.0};
        double chargeTimeSumC{0.0};
    };

    std::vector<Bucket> buckets;
    std::unordered_map<std::uint64_t, std::vector<std::size_t>> bucketsByCell;

    for (const auto& c : contributions) {
        if (!(c.inducedChargeC > 0.0)) {
            continue;
        }

        auto& indices = bucketsByCell[c.cellID];
        Bucket* bucket = nullptr;
        for (std::size_t idx : indices) {
            if (std::abs(buckets[idx].anchorTimeNs - c.timeNs) < integrationWindowNs) {
                bucket = &buckets[idx];
                break;
            }
        }
        if (bucket == nullptr) {
            Bucket fresh{};
            fresh.channel.cellID = c.cellID;
            fresh.channel.position = c.position;
            fresh.anchorTimeNs = c.timeNs;
            indices.push_back(buckets.size());
            buckets.push_back(fresh);
            bucket = &buckets.back();
        }

        bucket->channel.chargeC += c.inducedChargeC;
        bucket->channel.energyGeV += c.inducedEnergyGeV;
        bucket->chargeTimeSumC += c.inducedChargeC * c.timeNs;
        bucket->channel.contributors.emplace_back(c.simHitIndex, c.inducedChargeC);
    }

    std::vector<AggregatedChannel> channels;
    channels.reserve(buckets.size());
    for (auto& bucket : buckets) {
        AggregatedChannel channel = std::move(bucket.channel);
        channel.timeNs =
            (channel.chargeC > 0.0) ? bucket.chargeTimeSumC / channel.chargeC : bucket.anchorTimeNs;
        // Convert stored induced charge into a normalized truth weight.
        for (auto& [simHitIndex, weight] : channel.contributors) {
            weight = (channel.chargeC > 0.0) ? weight / channel.chargeC : 0.0;
        }
        channels.push_back(std::move(channel));
    }

    return channels;
}

// ---------------------------------------------------------------------------
// Per-hit processing (public for tests).
// ---------------------------------------------------------------------------
LGADChargeSharingRecon::SingleHitResult
LGADChargeSharingRecon::processSingleHit(const SingleHitInput& input) const {
    return processSingleHitImpl(input);
}

LGADChargeSharingRecon::SingleHitResult
LGADChargeSharingRecon::processSingleHitImpl(const SingleHitInput& input) const {
    SingleHitResult result{};

    PixelLocation nearest{};
    if (input.pixelIndexHint.has_value() && input.pixelHintMM.has_value()) {
        nearest.center = *input.pixelHintMM;
        nearest.indexI = input.pixelIndexHint->first;
        nearest.indexJ = input.pixelIndexHint->second;
    } else if (input.pixelIndexHint.has_value()) {
        const auto& [idxI, idxJ] = *input.pixelIndexHint;
        nearest = pixelLocationFromIndices(idxI, idxJ);
    } else if (input.pixelHintMM.has_value()) {
        nearest = findNearestPixelFallback(*input.pixelHintMM);
    } else {
        nearest = findNearestPixelFallback(input.hitPositionMM);
    }

    result.nearestPixelCenterMM = nearest.center;
    result.pixelRowIndex = nearest.indexI;
    result.pixelColIndex = nearest.indexJ;

    const double hitX = input.hitPositionMM[0];
    const double hitY = input.hitPositionMM[1];

    core::NeighborhoodConfig neighborCfg;
    neighborCfg.activeMode = core::ActivePixelMode::Neighborhood;
    neighborCfg.radius = m_cfg.neighborhoodRadius;
    neighborCfg.pixelSizeMM = m_geom.pixelSizeXMM;
    neighborCfg.pixelSizeYMM = m_geom.pixelSizeYMM;
    neighborCfg.pixelSpacingMM = m_geom.pixelSpacingXMM;
    neighborCfg.pixelSpacingYMM = m_geom.pixelSpacingYMM;
    neighborCfg.d0Micron = m_cfg.d0Micron;
    if (input.indexBounds.has_value()) {
        const auto& bounds = *input.indexBounds;
        neighborCfg.minIndexX = bounds.minX;
        neighborCfg.minIndexY = bounds.minY;
        neighborCfg.numPixelsX = bounds.maxX - bounds.minX + 1;
        neighborCfg.numPixelsY = bounds.maxY - bounds.minY + 1;
    } else {
        neighborCfg.numPixelsX = m_geom.numPixelsX;
        neighborCfg.numPixelsY = m_geom.numPixelsY;
        neighborCfg.minIndexX = m_geom.hasBoundsX() ? m_geom.minIndexX : 0;
        neighborCfg.minIndexY = m_geom.hasBoundsY() ? m_geom.minIndexY : 0;
    }

    core::NeighborhoodResult neighborhood = core::calculateNeighborhood(
        hitX, hitY, nearest.indexI, nearest.indexJ, nearest.center[0], nearest.center[1], neighborCfg);

    const double edepEV = input.energyDepositGeV * 1.0e9;
    const double numElectrons = (m_cfg.ionizationEnergyEV > 0.0) ? (edepEV / m_cfg.ionizationEnergyEV) : 0.0;
    const double totalChargeElectrons = numElectrons * m_cfg.amplificationFactor;
    const double totalChargeCoulombs = totalChargeElectrons * kElementaryChargeC;

    // Induced (noiseless) charge per pad. Electronic noise is applied later,
    // once per aggregated channel, in process().
    for (auto& pixel : neighborhood.pixels) {
        if (pixel.inBounds) {
            pixel.charge = pixel.fraction * totalChargeCoulombs;
        }
    }

    result.neighbors.reserve(neighborhood.pixels.size());
    for (const auto& pixel : neighborhood.pixels) {
        if (!pixel.inBounds)
            continue;
        NeighborData n{};
        n.fraction = pixel.fraction;
        n.chargeC = pixel.charge;
        n.distanceMM = pixel.distance;
        n.alphaRad = pixel.alpha;
        n.pixelXMM = pixel.centerX;
        n.pixelYMM = pixel.centerY;
        n.di = pixel.di;
        n.dj = pixel.dj;
        result.neighbors.push_back(n);
    }

    result.totalCollectedChargeC = totalChargeCoulombs;
    result.neighborhoodRadius = m_cfg.neighborhoodRadius;
    result.numActiveNeighbors = 0;
    for (const auto& n : result.neighbors) {
        if (n.chargeC > 0.0) {
            ++result.numActiveNeighbors;
        }
    }

    return result;
}

// ---------------------------------------------------------------------------
// Pixel location helpers
// ---------------------------------------------------------------------------
LGADChargeSharingRecon::PixelLocation
LGADChargeSharingRecon::findNearestPixelFallback(const std::array<double, 3>& positionMM) const {
    int i = positionToIndex(positionMM[0], m_geom.pixelSpacingXMM, m_geom.gridOffsetXMM);
    int j = positionToIndex(positionMM[1], m_geom.pixelSpacingYMM, m_geom.gridOffsetYMM);

    if (m_geom.hasBoundsX()) {
        i = std::clamp(i, m_geom.minIndexX, m_geom.maxIndexX);
    }
    if (m_geom.hasBoundsY()) {
        j = std::clamp(j, m_geom.minIndexY, m_geom.maxIndexY);
    }
    return pixelLocationFromIndices(i, j);
}

LGADChargeSharingRecon::PixelLocation
LGADChargeSharingRecon::pixelLocationFromIndices(int indexI, int indexJ) const {
    PixelLocation loc{};
    loc.indexI = indexI;
    loc.indexJ = indexJ;
    const double pixelCenterZ =
        m_geom.detectorZCenterMM + m_geom.detectorThicknessMM / 2.0 + m_geom.pixelThicknessMM / 2.0;
    loc.center = {indexToPosition(indexI, m_geom.pixelSpacingXMM, m_geom.gridOffsetXMM),
                  indexToPosition(indexJ, m_geom.pixelSpacingYMM, m_geom.gridOffsetYMM), pixelCenterZ};
    return loc;
}

} // namespace eicrecon
