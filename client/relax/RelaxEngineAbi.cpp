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
#include "eon/relax/eon_relax_engine.h"

#include "relax/HostSurfacePotential.h"
#include "eon/Eigen.h"
#include "eon/Matter.h"
#include "eon/MinModeSaddleSearch.h"
#include "eon/NudgedElasticBand.h"
#include "eon/Parameters.h"
#include "version.h"

#include <cstdio>
#include <cstring>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>

namespace {

uint64_t fnv1a64(const char *s) {
  uint64_t h = 14695981039346656037ULL;
  if (!s) {
    return h;
  }
  for (const unsigned char *p = reinterpret_cast<const unsigned char *>(s); *p;
       ++p) {
    h ^= static_cast<uint64_t>(*p);
    h *= 1099511628211ULL;
  }
  return h;
}

void set_err(char *errbuf, size_t errlen, const char *msg) {
  if (!errbuf || errlen == 0) {
    return;
  }
  std::snprintf(errbuf, errlen, "%s", msg ? msg : "");
}

bool known_kind(eon_relax_kind_t kind) {
  return kind == EON_RELAX_KIND_NEB || kind == EON_RELAX_KIND_SADDLE;
}

eon_relax_status_t map_neb(NudgedElasticBand::NEBStatus st) {
  using S = NudgedElasticBand::NEBStatus;
  switch (st) {
  case S::GOOD:
    return EON_RELAX_CONVERGED;
  case S::BAD_MAX_ITERATIONS:
    return EON_RELAX_MAX_ITERATIONS;
  case S::MAX_UNCERTAINTY:
    return EON_RELAX_MAX_UNCERTAINTY;
  case S::INIT:
    return EON_RELAX_INIT;
  case S::RUNNING:
    return EON_RELAX_OK;
  }
  return EON_RELAX_INVALID_PARAMETER;
}

void fill_matter(Matter &m, const eon_relax_band_t *band, long image,
                 uint64_t epoch) {
  m.resize(band->n_atoms);
  Eigen::VectorXi z(band->n_atoms);
  for (long a = 0; a < band->n_atoms; ++a) {
    z(a) = band->atomic_nrs[a];
  }
  m.setAtomicNrs(z);
  m.setMasses(Eigen::VectorXd::Ones(band->n_atoms));
  AtomMatrix pos(band->n_atoms, 3);
  const double *src = band->positions + image * 3 * band->n_atoms;
  // Match Potential::force / Matter storage (Eigen column-major n x 3).
  Eigen::Map<const AtomMatrix> mapped(src, band->n_atoms, 3);
  pos = mapped;
  m.setPositions(pos);
  Matrix3d cell;
  const double *b = band->boxes + image * 9;
  cell << b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8];
  m.setCell(cell);
  if (band->is_fixed) {
    for (long a = 0; a < band->n_atoms; ++a) {
      m.setFixed(a, band->is_fixed[a] ? 1 : 0);
    }
  }
  m.setSurfaceEpoch(epoch);
}

void write_image(double *dst, const Matter &m, long n_atoms) {
  const AtomMatrix &pos = m.getPositions();
  Eigen::Map<AtomMatrix> mapped(dst, n_atoms, 3);
  mapped = pos;
}

} // namespace

struct EonRelaxEngine {
  Parameters params;
  eon_relax_kind_t kind{EON_RELAX_KIND_NEB};
  uint64_t epoch{0};
};

extern "C" {

int eon_relax_abi_version(void) { return EON_RELAX_ABI_VERSION; }

int eon_relax_abi_stamp(eon_relax_abi_stamp_t *out) {
  if (!out) {
    return EON_RELAX_INVALID_PARAMETER;
  }
  out->major = EON_RELAX_ABI_MAJOR;
  out->minor = EON_RELAX_ABI_MINOR;
  out->layout_revision = EON_RELAX_ABI_LAYOUT;
  return EON_RELAX_OK;
}

uint64_t eon_relax_version_hash(void) {
  const uint64_t abi = (uint64_t)EON_RELAX_ABI_MAJOR << 48 |
                       (uint64_t)EON_RELAX_ABI_MINOR << 32 |
                       (uint64_t)EON_RELAX_ABI_LAYOUT << 16;
  return abi ^ fnv1a64(GIT_HASH.c_str());
}

EonRelaxEngine *eon_relax_create(const void *config, size_t config_len,
                                 char *errbuf, size_t errlen) {
  if (config != nullptr && config_len > 0) {
    set_err(errbuf, errlen,
            "RelaxEngineParams capnp parse is not wired; pass NULL config");
    return nullptr;
  }
  auto *eng = new (std::nothrow) EonRelaxEngine();
  if (!eng) {
    set_err(errbuf, errlen, "out of memory");
    return nullptr;
  }
  // Host owns endpoint relaxation; keep NSDMI otherwise.
  eng->params.neb_options.endpoints.minimize = false;
  return eng;
}

int eon_relax_set_kind(EonRelaxEngine *eng, eon_relax_kind_t kind) {
  if (!eng) {
    return EON_RELAX_INVALID_PARAMETER;
  }
  if (!known_kind(kind)) {
    return EON_RELAX_UNKNOWN_KIND;
  }
  eng->kind = kind;
  return EON_RELAX_OK;
}

int eon_relax_set_surface_epoch(EonRelaxEngine *eng, uint64_t epoch) {
  if (!eng) {
    return EON_RELAX_INVALID_PARAMETER;
  }
  eng->epoch = epoch;
  return EON_RELAX_OK;
}

int eon_relax_run(EonRelaxEngine *eng, const eon_relax_band_t *band,
                  eon_relax_surface_fn surface, void *user,
                  eon_relax_outcome_t *out) {
  if (!eng || !band || !surface || !out) {
    return EON_RELAX_INVALID_PARAMETER;
  }
  if (band->n_atoms <= 0 || !band->positions || !band->atomic_nrs ||
      !band->boxes) {
    return EON_RELAX_INVALID_PARAMETER;
  }
  if (!known_kind(eng->kind)) {
    return EON_RELAX_UNKNOWN_KIND;
  }
  std::memset(out, 0, sizeof(*out));
  out->version_hash = eon_relax_version_hash();
  out->surface_epoch = eng->epoch;

  auto pot = std::make_shared<eonc::HostSurfacePotential>(surface, user,
                                                          eng->epoch);
  try {
    if (eng->kind == EON_RELAX_KIND_NEB) {
      if (band->n_images < 2) {
        return EON_RELAX_INVALID_PARAMETER;
      }
      auto initial = std::make_shared<Matter>(pot, eng->params);
      auto final = std::make_shared<Matter>(pot, eng->params);
      fill_matter(*initial, band, 0, eng->epoch);
      fill_matter(*final, band, band->n_images - 1, eng->epoch);
      std::unique_ptr<NudgedElasticBand> neb;
      if (band->n_images == 2) {
        neb = std::make_unique<NudgedElasticBand>(initial, final, eng->params,
                                                  pot);
      } else {
        eng->params.neb_options.image_count = band->n_images - 2;
        std::vector<Matter> path;
        path.reserve(static_cast<size_t>(band->n_images));
        for (long i = 0; i < band->n_images; ++i) {
          Matter img(pot, eng->params);
          fill_matter(img, band, i, eng->epoch);
          path.push_back(std::move(img));
        }
        neb = std::make_unique<NudgedElasticBand>(std::move(path), eng->params,
                                                  pot);
      }
      const auto st = neb->compute();
      out->status = map_neb(st);
      out->climbing_image = neb->climbingImage;
      out->max_force = neb->convergenceForce();
      out->surface_epoch = pot->surfaceEpoch();
      if (band->n_images == 2) {
        write_image(band->positions, *neb->path.front(), band->n_atoms);
        write_image(band->positions + 3 * band->n_atoms, *neb->path.back(),
                    band->n_atoms);
      } else {
        for (long i = 0; i < band->n_images; ++i) {
          write_image(band->positions + i * 3 * band->n_atoms, *neb->path[i],
                      band->n_atoms);
        }
      }
      return out->status >= 0 ? EON_RELAX_OK : out->status;
    }

    if (band->n_images < 1 || !band->mode) {
      return EON_RELAX_INVALID_PARAMETER;
    }
    auto matter = std::make_shared<Matter>(pot, eng->params);
    fill_matter(*matter, band, 0, eng->epoch);
    AtomMatrix mode(band->n_atoms, 3);
    Eigen::Map<const AtomMatrix> mmap(band->mode, band->n_atoms, 3);
    mode = mmap;
    const double reactant_e = matter->getPotentialEnergy();
    MinModeSaddleSearch search(matter, mode, reactant_e, eng->params, pot);
    const int sst = search.run();
    write_image(band->positions, *matter, band->n_atoms);
    out->iterations = search.getIterationCount();
    out->surface_epoch = pot->surfaceEpoch();
    if (sst == MinModeSaddleSearch::STATUS_GOOD) {
      out->status = EON_RELAX_CONVERGED;
      return EON_RELAX_OK;
    }
    if (sst == MinModeSaddleSearch::STATUS_BAD_MAX_ITERATIONS) {
      out->status = EON_RELAX_MAX_ITERATIONS;
      return EON_RELAX_OK;
    }
    out->status = EON_RELAX_UNAVAILABLE;
    return EON_RELAX_OK;
  } catch (const std::exception &) {
    out->status = EON_RELAX_SURFACE_FAILED;
    return EON_RELAX_SURFACE_FAILED;
  }
}

void eon_relax_destroy(EonRelaxEngine *eng) { delete eng; }

} // extern "C"
