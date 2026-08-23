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
#include "TestUtils.hpp"
#include "catch2/catch_amalgamated.hpp"
#include "eon/ConFileIO.h"
#include "eon/DynLib.h"
#include "eon/Eigen.h"
#include "eon/HelperFunctions.h"
#include "eon/Matter.h"
#include "eon/Parameters.h"
#include "eon/Potential.h"
#include "eon/relax/eon_relax_engine.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace tests {

static eonc::helpers::test::QuillTestLogger _quill_setup;

struct SurfaceCtx {
  std::shared_ptr<Potential> pot;
  double variance{0.0};
  uint64_t epoch{0};
  long calls{0};
};

static int surface_forward(void *user, long n_images, long n_atoms,
                           const double *positions, const int *atomic_nrs,
                           const double *boxes, const long *,
                           double *energies, double *forces, double *variances,
                           uint64_t *epoch_out) {
  auto *ctx = static_cast<SurfaceCtx *>(user);
  ++ctx->calls;
  for (long s = 0; s < n_images; ++s) {
    double var = 0.0;
    ctx->pot->force(n_atoms, positions + s * 3 * n_atoms, atomic_nrs,
                    forces + s * 3 * n_atoms, &energies[s], &var,
                    boxes + s * 9);
    if (variances) {
      variances[s] = ctx->variance;
    }
  }
  if (epoch_out) {
    *epoch_out = ctx->epoch;
  }
  return 0;
}

static int surface_fail(void *, long, long, const double *, const int *,
                        const double *, const long *, double *, double *,
                        double *, uint64_t *) {
  return -3;
}

static void pack_image(std::vector<double> &pos, std::vector<double> &boxes,
                       std::vector<int> &z, const Matter &m, long image) {
  const long n = m.numberOfAtoms();
  const AtomMatrix &p = m.getPositions();
  double *dst = pos.data() + image * 3 * n;
  Eigen::Map<AtomMatrix> mapped(dst, n, 3);
  mapped = p;
  const Matrix3d cell = m.getCell();
  double *b = boxes.data() + image * 9;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      b[i * 3 + j] = cell(i, j);
    }
  }
  if (image == 0) {
    auto nrs = m.getAtomicNrs();
    z.resize(static_cast<size_t>(n));
    for (long a = 0; a < n; ++a) {
      z[static_cast<size_t>(a)] = nrs(a);
    }
  }
}

static void pack_linear_band(std::vector<double> &pos, std::vector<double> &boxes,
                             std::vector<int> &z, const Matter &a,
                             const Matter &b, long n_images) {
  const long n = a.numberOfAtoms();
  pos.assign(static_cast<size_t>(n_images * 3 * n), 0.0);
  boxes.assign(static_cast<size_t>(n_images * 9), 0.0);
  pack_image(pos, boxes, z, a, 0);
  pack_image(pos, boxes, z, b, n_images - 1);
  const AtomMatrix pa = a.getPositions();
  const AtomMatrix pb = b.getPositions();
  for (long im = 1; im < n_images - 1; ++im) {
    const double t = static_cast<double>(im) / static_cast<double>(n_images - 1);
    AtomMatrix mid = pa + t * (pb - pa);
    double *dst = pos.data() + im * 3 * n;
    Eigen::Map<AtomMatrix> mapped(dst, n, 3);
    mapped = mid;
    std::memcpy(boxes.data() + im * 9, boxes.data(), 9 * sizeof(double));
  }
}

TEST_CASE("relax engine ABI stamp is layout 1.0.1", "[relax][abi]") {
  REQUIRE(eon_relax_abi_version() == static_cast<int>(EON_RELAX_ABI_VERSION));
  REQUIRE(EON_RELAX_ABI_UNPACK_MAJOR(eon_relax_abi_version()) == 1);
  REQUIRE(EON_RELAX_ABI_UNPACK_MINOR(eon_relax_abi_version()) == 0);
  REQUIRE(EON_RELAX_ABI_UNPACK_LAYOUT(eon_relax_abi_version()) == 1);
  REQUIRE(eon_relax_available() == 1);
  eon_relax_abi_stamp_t stamp{};
  REQUIRE(eon_relax_abi_stamp(&stamp) == EON_RELAX_OK);
  REQUIRE(stamp.major == 1);
  REQUIRE(stamp.minor == 0);
  REQUIRE(stamp.layout_revision == 1);
  REQUIRE(eon_relax_abi_stamp(nullptr) == EON_RELAX_INVALID_PARAMETER);
  REQUIRE(eon_relax_version_hash() != 0);
}

TEST_CASE("relax engine destroy NULL is a no-op", "[relax][abi]") {
  eon_relax_destroy(nullptr);
}

TEST_CASE("relax engine status_name is fail-closed", "[relax][abi]") {
  REQUIRE(std::string(eon_relax_status_name(EON_RELAX_KIND_NEB,
                                            EON_RELAX_NEB_GOOD)) == "GOOD");
  REQUIRE(std::string(eon_relax_status_name(
              EON_RELAX_KIND_NEB, EON_RELAX_NEB_MAX_UNCERTAINTY)) ==
          "MAX_UNCERTAINTY");
  REQUIRE(eon_relax_status_name(EON_RELAX_KIND_NEB, 99) == nullptr);
  REQUIRE(eon_relax_status_name(EON_RELAX_KIND_INVALID, 0) == nullptr);
}

TEST_CASE("relax engine create NULL config and reject unknown capnp",
          "[relax][abi]") {
  char err[128]{};
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, err, sizeof(err));
  REQUIRE(eng != nullptr);
  const char blob[] = "not-capnp";
  REQUIRE(eon_relax_create(blob, sizeof(blob), err, sizeof(err)) == nullptr);
  REQUIRE(std::string(err).find("-13") != std::string::npos);
  REQUIRE(eon_relax_run(nullptr, nullptr, nullptr, nullptr, nullptr) ==
          EON_RELAX_NULL_ENGINE);
  eon_relax_band_t short_band{};
  short_band.n_images = 2;
  short_band.n_atoms = 1;
  double pos1[6]{};
  int z1[1]{1};
  double box1[18]{};
  short_band.positions = pos1;
  short_band.atomic_nrs = z1;
  short_band.boxes = box1;
  eon_relax_outcome_t out{};
  REQUIRE(eon_relax_run(eng, &short_band, surface_fail, nullptr, &out) ==
          EON_RELAX_BAND_TOO_SHORT);
  short_band.n_images = 4;
  REQUIRE(eon_relax_run(eng, &short_band, surface_fail, nullptr, &out) ==
          EON_RELAX_BAND_SIZE);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine NEB on LJ cluster reaches a finite force",
          "[relax][neb]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto reactant = std::make_shared<Matter>(pot, params);
  auto product = std::make_shared<Matter>(pot, params);
  REQUIRE(reactant->con2matter(std::string("reactant.con")) ==
          eonc::io::IoStatus::Ok);
  REQUIRE(product->con2matter(std::string("reactant.con")) ==
          eonc::io::IoStatus::Ok);
  auto p = product->getPositions();
  p(0, 0) += 0.5;
  p(0, 1) -= 0.3;
  p(0, 2) += 0.2;
  product->setPositions(p);

  const long n = reactant->numberOfAtoms();
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int> z;
  pack_linear_band(pos, boxes, z, *reactant, *product, n_images);

  SurfaceCtx ctx{pot, 0.0, 0, 0};
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.n_images = n_images;
  band.n_atoms = n;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  eon_relax_outcome_t out{};
  const int rc = eon_relax_run(eng, &band, surface_forward, &ctx, &out);
  REQUIRE(rc == EON_RELAX_OK);
  REQUIRE(out.status == EON_RELAX_NEB_GOOD);
  REQUIRE(out.kind == EON_RELAX_KIND_NEB);
  REQUIRE(std::isfinite(out.max_force));
  REQUIRE(out.max_force <= 0.01 + 1e-6);
  REQUIRE(ctx.calls > 0);
  REQUIRE(out.version_hash != 0);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine MAX_UNCERTAINTY surfaces from host variance",
          "[relax][neb][uncertainty]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto reactant = std::make_shared<Matter>(pot, params);
  auto product = std::make_shared<Matter>(pot, params);
  REQUIRE(reactant->con2matter(std::string("reactant.con")) ==
          eonc::io::IoStatus::Ok);
  REQUIRE(product->con2matter(std::string("reactant.con")) ==
          eonc::io::IoStatus::Ok);
  auto p = product->getPositions();
  p(0, 0) += 0.5;
  product->setPositions(p);

  const long n = reactant->numberOfAtoms();
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int> z;
  pack_linear_band(pos, boxes, z, *reactant, *product, n_images);

  SurfaceCtx ctx{pot, 10.0, 0, 0};
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.n_images = n_images;
  band.n_atoms = n;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  eon_relax_outcome_t out{};
  const int rc = eon_relax_run(eng, &band, surface_forward, &ctx, &out);
  REQUIRE(rc == EON_RELAX_OK);
  REQUIRE(out.status == EON_RELAX_NEB_MAX_UNCERTAINTY);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine surface failure is fail-closed", "[relax][abi]") {
  Parameters params;
  params.potential_options.potential = PotType::LJ;
  auto pot = eonc::helpers::makePotential(PotType::LJ, params);
  auto reactant = std::make_shared<Matter>(pot, params);
  REQUIRE(reactant->con2matter(std::string("reactant.con")) ==
          eonc::io::IoStatus::Ok);
  const long n = reactant->numberOfAtoms();
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int> z;
  pack_linear_band(pos, boxes, z, *reactant, *reactant, n_images);

  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.n_images = n_images;
  band.n_atoms = n;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();
  eon_relax_outcome_t out{};
  const int rc = eon_relax_run(eng, &band, surface_fail, nullptr, &out);
  REQUIRE(rc == EON_RELAX_SURFACE_FATAL);
  eon_relax_destroy(eng);
}

TEST_CASE("libeon_relax_engine exports the C waist", "[relax][dlopen]") {
  const char *env = std::getenv("EON_RELAX_ENGINE");
  std::vector<std::string> candidates;
  if (env && *env) {
    candidates.emplace_back(env);
  }
  candidates.emplace_back("libeon_relax_engine.so");
  eonc::dynlib::Handle h{};
  for (const auto &path : candidates) {
    h = eonc::dynlib::open(path.c_str());
    if (h) {
      break;
    }
  }
  REQUIRE(h);
  REQUIRE(eonc::dynlib::sym(h, "eon_relax_create") != nullptr);
  REQUIRE(eonc::dynlib::sym(h, "eon_relax_run") != nullptr);
  REQUIRE(eonc::dynlib::sym(h, "eon_relax_destroy") != nullptr);
  REQUIRE(eonc::dynlib::sym(h, "eon_relax_abi_version") != nullptr);
  auto abi = eonc::dynlib::loadSym<int (*)()>(h, "eon_relax_abi_version");
  REQUIRE(abi != nullptr);
  REQUIRE(abi() == EON_RELAX_ABI_VERSION);
  eonc::dynlib::close(h);
}

} // namespace tests
