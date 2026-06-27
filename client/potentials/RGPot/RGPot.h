/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#pragma once

#include "../../Potential.h"
#include <memory>
#include <string>

class RGPotEngine;

/**
 * Potential backed by rgpot NWChemPot / CPMDPot (dlopen C ABI engines).
 * Configure via [RGPot] INI; energy eV, forces eV/Angstrom.
 */
class RGPot final : public Potential {
public:
  explicit RGPot(const Parameters &p);
  ~RGPot() override;

  RGPot(const RGPot &) = delete;
  RGPot &operator=(const RGPot &) = delete;

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }
  [[nodiscard]] const std::string &backend() const noexcept { return backend_; }
  [[nodiscard]] bool engineAvailable() const;

private:
  std::unique_ptr<RGPotEngine> impl_;
  std::string backend_;
};
