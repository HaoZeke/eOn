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
 * Potential backed by rgpot NWChemPot / CPMDPot (in-process dlopen of
 * libnwchemc / libcpmdc). Not potserv RPC. Configure via [RgpotPot] INI;
 * energy eV, forces eV/Angstrom.
 */
class RgpotPot final : public Potential {
public:
  explicit RgpotPot(const Parameters &p);
  ~RgpotPot() override;

  RgpotPot(const RgpotPot &) = delete;
  RgpotPot &operator=(const RgpotPot &) = delete;

  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

  [[nodiscard]] bool isThreadSafe() const noexcept override { return false; }
  /// NWChem molecular SCF does not support PBC; CPMD may be periodic.
  [[nodiscard]] bool requiresIsolatedMoleculeLayout() const noexcept override {
    return backend_ == "nwchemc";
  }
  [[nodiscard]] const std::string &backend() const noexcept { return backend_; }
  [[nodiscard]] bool engineAvailable() const;

  /// Block on multi-rank force collectives (non-zero MPI ranks only).
  void mpiWorkerLoop() const;

  static int mpiWorldSize();
  static int mpiWorldRank();

  /**
   * Multi-rank driver pot (rank 0). Kept alive for the whole job so we do not
   * destroy/finalize libnwchemc while workers sit in mpiWorkerLoop (that is an
   * MPI deadlock: finalize/GA collectives vs Bcast of N).
   * Potential factory returns this instance instead of constructing a second
   * pot on rank 0.
   */
  static void setMpiDriverInstance(std::shared_ptr<RgpotPot> pot);
  static std::shared_ptr<RgpotPot> mpiDriverInstance();
  static void clearMpiDriverInstance();

private:
  std::unique_ptr<RGPotEngine> impl_;
  std::string backend_;
};
