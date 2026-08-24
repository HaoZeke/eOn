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

static int surface_forward(void *user, eon_relax_surface_request_t *req) {
  auto *ctx = static_cast<SurfaceCtx *>(user);
  ++ctx->calls;
  const long n_images = static_cast<long>(req->n_images);
  const long n_atoms = static_cast<long>(req->n_atoms);
  const double *positions = req->positions;
  const int32_t *atomic_nrs_i32 = req->atomic_nrs;
  std::vector<int> nrs(atomic_nrs_i32, atomic_nrs_i32 + n_atoms);
  const int *atomic_nrs = nrs.data();
  const double *boxes = req->boxes;
  double *energies = req->energies;
  double *forces = req->forces;
  double *variances = req->variances;
  uint64_t *epoch_out = req->epoch_out;
  for (long s = 0; s < n_images; ++s) {
    if (ctx->pot) {
      double var = 0.0;
      ctx->pot->force(n_atoms, positions + s * 3 * n_atoms, atomic_nrs,
                      forces + s * 3 * n_atoms, &energies[s], &var,
                      boxes + s * 9);
    } else {
      const double *p = positions + s * 3 * n_atoms;
      double *f = forces + s * 3 * n_atoms;
      double e = 0.0;
      for (long a = 0; a < n_atoms; ++a) {
        const double x = p[3 * a];
        const double y = p[3 * a + 1];
        const double z = p[3 * a + 2];
        e += 0.5 * (x * x + y * y + z * z);
        f[3 * a] = -x;
        f[3 * a + 1] = -y;
        f[3 * a + 2] = -z;
      }
      energies[s] = e;
    }
    if (variances) {
      variances[s] = ctx->variance;
    }
  }
  if (epoch_out) {
    *epoch_out = ctx->epoch;
  }
  return 0;
}

static void pack_harmonic_band(std::vector<double> &pos,
                               std::vector<double> &boxes, std::vector<int32_t> &z,
                               long n_images) {
  const long n_atoms = 1;
  pos.assign(static_cast<size_t>(n_images * 3 * n_atoms), 0.0);
  boxes.assign(static_cast<size_t>(n_images * 9), 0.0);
  z.assign(1, 1);
  for (long im = 0; im < n_images; ++im) {
    const double t =
        static_cast<double>(im) / static_cast<double>(n_images - 1);
    pos[static_cast<size_t>(im * 3)] = 0.3 - 0.6 * t;
    double *b = boxes.data() + im * 9;
    b[0] = 10.0;
    b[4] = 10.0;
    b[8] = 10.0;
  }
}

static int surface_fail(void *, eon_relax_surface_request_t *) { return -3; }

TEST_CASE("relax engine ABI stamp is 2.0", "[relax][abi]") {
  REQUIRE(eon_relax_abi_version() ==
          static_cast<int>((EON_RELAX_ABI_MAJOR << 16) | EON_RELAX_ABI_MINOR));
  REQUIRE(eon_relax_available() == 1);
  eon_relax_version_t stamp{};
  REQUIRE(eon_relax_abi_stamp(&stamp) == EON_RELAX_OK);
  REQUIRE(stamp.major == 2);
  REQUIRE(stamp.minor == 0);
  REQUIRE(eon_relax_abi_stamp(nullptr) == EON_RELAX_INVALID_PARAMETER);
  REQUIRE(eon_relax_version_hash() != 0);
  const char *id = eon_relax_version_hash_str();
  REQUIRE(id != nullptr);
  REQUIRE(std::string(id).find("+git.") != std::string::npos);
}

TEST_CASE("relax engine destroy NULL is a no-op", "[relax][abi]") {
  eon_relax_destroy(nullptr);
  REQUIRE(eon_relax_last_error(nullptr) == EON_RELAX_NULL_ENGINE);
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
  short_band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  short_band.n_images = 2;
  short_band.n_atoms = 1;
  double pos1[6]{};
  int32_t z1[1]{1};
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
  REQUIRE(eon_relax_last_error(eng) == EON_RELAX_BAND_SIZE);
  short_band.n_images = 7;
  short_band.n_atoms = 0;
  REQUIRE(eon_relax_run(eng, &short_band, surface_fail, nullptr, &out) ==
          EON_RELAX_NATOMS);
  short_band.n_atoms = 1;
  short_band.positions = nullptr;
  REQUIRE(eon_relax_run(eng, &short_band, surface_fail, nullptr, &out) ==
          EON_RELAX_NULL_POSITIONS);
  short_band.positions = pos1;
  short_band.atomic_nrs = nullptr;
  REQUIRE(eon_relax_run(eng, &short_band, surface_fail, nullptr, &out) ==
          EON_RELAX_NULL_ATOMIC_NRS);
  short_band.atomic_nrs = z1;
  short_band.boxes = nullptr;
  REQUIRE(eon_relax_run(eng, &short_band, surface_fail, nullptr, &out) ==
          EON_RELAX_NULL_BOXES);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine NEB on a harmonic surface converges",
          "[relax][neb]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);

  SurfaceCtx ctx{nullptr, 0.0, 0, 0};
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
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
  // Endpoints stay put: NudgedElasticBandJob owns endpoint min, not compute().
  REQUIRE(std::abs(pos[0] - 0.3) < 1e-6);
  REQUIRE(std::abs(pos[3]) <= 0.2 + 1e-6);
  eon_relax_destroy(eng);
}

TEST_CASE("relax engine MAX_UNCERTAINTY surfaces from host variance",
          "[relax][neb][uncertainty]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);

  SurfaceCtx ctx{nullptr, 10.0, 0, 0};
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
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
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);

  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
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
  REQUIRE(eonc::dynlib::sym(h, "eon_relax_status_name") != nullptr);
  REQUIRE(eonc::dynlib::sym(h, "eon_relax_last_error") != nullptr);
  REQUIRE(eonc::dynlib::sym(h, "eon_relax_version_hash_str") != nullptr);
  REQUIRE(eonc::dynlib::sym(h, "eon_relax_abi_version") != nullptr);
  auto abi = eonc::dynlib::loadSym<int (*)()>(h, "eon_relax_abi_version");
  REQUIRE(abi != nullptr);
  REQUIRE(abi() == static_cast<int>((EON_RELAX_ABI_MAJOR << 16) |
                                    EON_RELAX_ABI_MINOR));
  eonc::dynlib::close(h);
}


TEST_CASE("relax engine stepper converges the harmonic band one step at a "
          "time",
          "[relax][neb][stepper]") {
  const long n_images = 7;
  std::vector<double> pos;
  std::vector<double> boxes;
  std::vector<int32_t> z;
  pack_harmonic_band(pos, boxes, z, n_images);

  SurfaceCtx ctx{nullptr, 0.0, 0, 0};
  EonRelaxEngine *eng = eon_relax_create(nullptr, 0, nullptr, 0);
  REQUIRE(eng != nullptr);
  eon_relax_band_t band{};
  band.version = eon_relax_version_t EON_RELAX_VERSION_INIT;
  band.n_images = n_images;
  band.n_atoms = 1;
  band.positions = pos.data();
  band.atomic_nrs = z.data();
  band.boxes = boxes.data();

  eon_relax_outcome_t out{};
  int steps = 0;
  int rc = EON_RELAX_OK;
  do {
    rc = eon_relax_step(eng, &band, surface_forward, &ctx, &out);
    REQUIRE(rc == EON_RELAX_OK);
    ++steps;
  } while (out.status == EON_RELAX_NEB_RUNNING && steps < 500);
  REQUIRE(out.status == EON_RELAX_NEB_GOOD);
  REQUIRE(out.iterations <= steps);
  REQUIRE(std::isfinite(out.max_force));
  // Endpoints are the caller's; the stepper never moves them.
  REQUIRE(std::abs(pos[0] - 0.3) < 1e-12);
  REQUIRE(std::abs(pos[(n_images - 1) * 3] + 0.3) < 1e-12);

  // A different surface pointer without a reset is refused.
  REQUIRE(eon_relax_step(eng, &band, surface_fail, &ctx, &out) ==
          EON_RELAX_INVALID_PARAMETER);
  // Reset drops the band state; the next step reinitializes.
  REQUIRE(eon_relax_reset(eng) == EON_RELAX_OK);
  REQUIRE(eon_relax_step(eng, &band, surface_forward, &ctx, &out) ==
          EON_RELAX_OK);
  eon_relax_destroy(eng);
}

} // namespace tests
