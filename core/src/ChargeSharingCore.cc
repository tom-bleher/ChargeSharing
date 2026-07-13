// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024-2026 Tom Bleher, Igor Korover

#include "chargesharing/core/ChargeSharingCore.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace chargesharing::core {

std::vector<const NeighborPixel*> NeighborhoodResult::getCenterRow() const {
    std::vector<const NeighborPixel*> row;
    for (const auto& p : pixels) {
        if (p.dj == 0 && p.inBounds) {
            row.push_back(&p);
        }
    }
    std::sort(row.begin(), row.end(),
              [](const NeighborPixel* a, const NeighborPixel* b) { return a->di < b->di; });
    return row;
}

std::vector<const NeighborPixel*> NeighborhoodResult::getCenterCol() const {
    std::vector<const NeighborPixel*> col;
    for (const auto& p : pixels) {
        if (p.di == 0 && p.inBounds) {
            col.push_back(&p);
        }
    }
    std::sort(col.begin(), col.end(),
              [](const NeighborPixel* a, const NeighborPixel* b) { return a->dj < b->dj; });
    return col;
}

NeighborhoodResult calculateNeighborhood(double hitX, double hitY, int centerI, int centerJ, double centerX,
                                         double centerY, const NeighborhoodConfig& config) {

    NeighborhoodResult result;
    result.centerPixelI = centerI;
    result.centerPixelJ = centerJ;
    result.centerPixelX = centerX;
    result.centerPixelY = centerY;

    const int radius = std::max(0, config.radius);
    const int gridDim = 2 * radius + 1;
    const double d0MM = config.d0Micron * constants::kMillimeterPerMicron;
    const double padW = config.pixelSizeMM;
    const double padH = (config.pixelSizeYMM > 0) ? config.pixelSizeYMM : padW;
    const double pitchX = config.pixelSpacingMM;
    const double pitchY = (config.pixelSpacingYMM > 0) ? config.pixelSpacingYMM : pitchX;

    result.pixels.reserve(gridDim * gridDim);

    // First pass: calculate weights for all pixels
    for (int di = -radius; di <= radius; ++di) {
        for (int dj = -radius; dj <= radius; ++dj) {
            NeighborPixel pixel;
            pixel.di = di;
            pixel.dj = dj;

            const int globalI = centerI + di;
            const int globalJ = centerJ + dj;

            // Check bounds if specified (indices can be negative for centered DD4hep grids)
            if (config.numPixelsX > 0 && (globalI < config.minIndexX || globalI >= config.minIndexX + config.numPixelsX)) {
                pixel.inBounds = false;
                result.pixels.push_back(pixel);
                continue;
            }
            if (config.numPixelsY > 0 && (globalJ < config.minIndexY || globalJ >= config.minIndexY + config.numPixelsY)) {
                pixel.inBounds = false;
                result.pixels.push_back(pixel);
                continue;
            }

            pixel.inBounds = true;
            pixel.centerX = centerX + di * pitchX;
            pixel.centerY = centerY + dj * pitchY;

            const double dx = hitX - pixel.centerX;
            const double dy = hitY - pixel.centerY;
            pixel.distance = calcDistanceToEdge(dx, dy, padW / 2.0, padH / 2.0);
            pixel.alpha = calcPadViewAngle(pixel.distance, padW, padH);

            pixel.weight = calcWeightLogA(pixel.distance, pixel.alpha, d0MM);

            result.pixels.push_back(pixel);
        }
    }

    // Determine which pixels contribute to denominator based on active mode
    std::vector<bool> isActive(result.pixels.size(), false);

    switch (config.activeMode) {
        case ActivePixelMode::Neighborhood:
            for (std::size_t i = 0; i < result.pixels.size(); ++i) {
                isActive[i] = result.pixels[i].inBounds;
            }
            break;

        case ActivePixelMode::RowCol:
            for (std::size_t i = 0; i < result.pixels.size(); ++i) {
                const auto& p = result.pixels[i];
                isActive[i] = p.inBounds && (p.di == 0 || p.dj == 0);
            }
            break;

        case ActivePixelMode::RowCol3x3:
            for (std::size_t i = 0; i < result.pixels.size(); ++i) {
                const auto& p = result.pixels[i];
                const bool inCross = (p.di == 0 || p.dj == 0);
                const bool in3x3 = (std::abs(p.di) <= 1 && std::abs(p.dj) <= 1);
                isActive[i] = p.inBounds && (inCross || in3x3);
            }
            break;

        case ActivePixelMode::ChargeBlock2x2:
        case ActivePixelMode::ChargeBlock3x3: {
            std::size_t maxIdx = 0;
            double maxWeight = -1.0;
            for (std::size_t i = 0; i < result.pixels.size(); ++i) {
                if (result.pixels[i].inBounds && result.pixels[i].weight > maxWeight) {
                    maxWeight = result.pixels[i].weight;
                    maxIdx = i;
                }
            }

            const int blockSize = (config.activeMode == ActivePixelMode::ChargeBlock3x3) ? 3 : 2;
            const int maxDi = result.pixels[maxIdx].di;
            const int maxDj = result.pixels[maxIdx].dj;

            int bestCornerDi = maxDi, bestCornerDj = maxDj;
            double bestSum = -1.0;
            const int maxOffset = (blockSize == 3) ? -2 : -1;

            for (int oi = maxOffset; oi <= 0; ++oi) {
                for (int oj = maxOffset; oj <= 0; ++oj) {
                    const int testCornerDi = maxDi + oi;
                    const int testCornerDj = maxDj + oj;

                    double sum = 0.0;
                    for (const auto& p : result.pixels) {
                        if (p.inBounds && p.di >= testCornerDi && p.di < testCornerDi + blockSize &&
                            p.dj >= testCornerDj && p.dj < testCornerDj + blockSize) {
                            sum += p.weight;
                        }
                    }

                    if (sum > bestSum) {
                        bestSum = sum;
                        bestCornerDi = testCornerDi;
                        bestCornerDj = testCornerDj;
                    }
                }
            }

            for (std::size_t i = 0; i < result.pixels.size(); ++i) {
                const auto& p = result.pixels[i];
                isActive[i] = p.inBounds && p.di >= bestCornerDi && p.di < bestCornerDi + blockSize &&
                              p.dj >= bestCornerDj && p.dj < bestCornerDj + blockSize;
            }
            break;
        }

        case ActivePixelMode::ThresholdAboveNoise:
            for (std::size_t i = 0; i < result.pixels.size(); ++i) {
                isActive[i] = result.pixels[i].inBounds;
            }
            break;
    }

    // Calculate total weight from active pixels
    result.totalWeight = 0.0;
    for (std::size_t i = 0; i < result.pixels.size(); ++i) {
        if (isActive[i]) {
            result.totalWeight += result.pixels[i].weight;
        }
    }

    // Normalize to get fractions (only active pixels get non-zero fractions)
    for (std::size_t i = 0; i < result.pixels.size(); ++i) {
        if (isActive[i] && result.totalWeight > 0.0) {
            result.pixels[i].fraction = result.pixels[i].weight / result.totalWeight;
        } else {
            result.pixels[i].fraction = 0.0;
        }
    }

    // -------------------------------------------------------------------------
    // Compute mode-specific fractions for diagnostics
    // These use different denominators: row-only, col-only, and block
    // -------------------------------------------------------------------------

    // Row-mode: denominator from center row only (dj=0)
    double rowWeightSum = 0.0;
    for (const auto& p : result.pixels) {
        if (p.inBounds && p.dj == 0) {
            rowWeightSum += p.weight;
        }
    }
    for (auto& p : result.pixels) {
        if (p.inBounds && rowWeightSum > 0.0) {
            p.fractionRow = p.weight / rowWeightSum;
        }
    }

    // Col-mode: denominator from center column only (di=0)
    double colWeightSum = 0.0;
    for (const auto& p : result.pixels) {
        if (p.inBounds && p.di == 0) {
            colWeightSum += p.weight;
        }
    }
    for (auto& p : result.pixels) {
        if (p.inBounds && colWeightSum > 0.0) {
            p.fractionCol = p.weight / colWeightSum;
        }
    }

    // Block-mode: find contiguous 2x2 block with highest weight sum
    std::size_t maxIdx = 0;
    double maxWeight = -1.0;
    for (std::size_t i = 0; i < result.pixels.size(); ++i) {
        if (result.pixels[i].inBounds && result.pixels[i].weight > maxWeight) {
            maxWeight = result.pixels[i].weight;
            maxIdx = i;
        }
    }

    const int blockSize = 2; // 2x2 block for diagnostics
    int blockCornerDi = result.pixels[maxIdx].di;
    int blockCornerDj = result.pixels[maxIdx].dj;
    double bestBlockSum = -1.0;

    for (int oi = -1; oi <= 0; ++oi) {
        for (int oj = -1; oj <= 0; ++oj) {
            const int testCornerDi = result.pixels[maxIdx].di + oi;
            const int testCornerDj = result.pixels[maxIdx].dj + oj;

            double sum = 0.0;
            for (const auto& p : result.pixels) {
                if (p.inBounds && p.di >= testCornerDi && p.di < testCornerDi + blockSize && p.dj >= testCornerDj &&
                    p.dj < testCornerDj + blockSize) {
                    sum += p.weight;
                }
            }

            if (sum > bestBlockSum) {
                bestBlockSum = sum;
                blockCornerDi = testCornerDi;
                blockCornerDj = testCornerDj;
            }
        }
    }

    for (auto& p : result.pixels) {
        if (p.inBounds && p.di >= blockCornerDi && p.di < blockCornerDi + blockSize && p.dj >= blockCornerDj &&
            p.dj < blockCornerDj + blockSize && bestBlockSum > 0.0) {
            p.fractionBlock = p.weight / bestBlockSum;
        }
    }

    return result;
}

} // namespace chargesharing::core
