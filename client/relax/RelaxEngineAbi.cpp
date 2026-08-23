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

bool known_kind(eon_relax_kind_t kind) { return EON_RELAX_KIND_IS_KNOWN(kind); }

int neb_status(NudgedElasticBand::NEBStatus st) {
  return static_cast<int>(st);
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

int eon_relax_abi_version(void) { return static_cast<int>(EON_RELAX_ABI_VERSION); }

int eon_relax_available(void) { return 1; }

int eon_relax_abi_stamp(eon_relax_abi_stamp_t *out) {
  if (!out) {
    return EON_RELAX_INVALID_PARAMETER;
  }
  out->major = EON_RELAX_ABI_MAJOR;
  out->minor = EON_RELAX_ABI_MINOR;
  out->layout_revision = EON_RELAX_ABI_LAYOUT;
  return EON_RELAX_OK;
}

const char *eon_relax_version_hash_str(void) {
  static const std::string id = VERSION + "+git." + GIT_HASH;
  return id.c_str();
}

uint64_t eon_relax_version_hash(void) {
  return fnv1a64(eon_relax_version_hash_str());
}

EonRelaxEngine *eon_relax_create(const void *config, size_t config_len,
                                 char *errbuf, size_t errlen) {
  if (config != nullptr && config_len > 0) {
    set_err(errbuf, errlen,
            "-13 RelaxEngineParams capnp parse is not wired; pass NULL config");
    return nullptr;
  }
  auto *eng = new (std::nothrow) EonRelaxEngine();
  if (!eng) {
    set_err(errbuf, errlen, "-19 out of memory");
    return nullptr;
  }
  eng->kind = EON_RELAX_KIND_NEB;
  return eng;
}

int eon_relax_set_surface_epoch(EonRelaxEngine *eng, uint64_t epoch) {
  if (!eng) {
    return EON_RELAX_INVALID_PARAMETER;
  }
  eng->epoch = epoch;
  return EON_RELAX_OK;
}

int eon_relax_run(EonRelaxEngine *eng, eon_relax_band_t *band,
                  eon_relax_surface_fn surface, void *user,
                  eon_relax_outcome_t *out) {
  if (!eng) {
    return EON_RELAX_NULL_ENGINE;
  }
  if (!band) {
    return EON_RELAX_NULL_BAND;
  }
  if (!surface) {
    return EON_RELAX_NULL_SURFACE;
  }
  if (!out) {
    return EON_RELAX_NULL_OUTCOME;
  }
  if (band->n_atoms <= 0) {
    return EON_RELAX_NATOMS;
  }
  if (!band->positions) {
    return EON_RELAX_NULL_POSITIONS;
  }
  if (!band->atomic_nrs) {
    return EON_RELAX_NULL_ATOMIC_NRS;
  }
  if (!band->boxes) {
    return EON_RELAX_NULL_BOXES;
  }
  if (!known_kind(eng->kind)) {
    return EON_RELAX_UNKNOWN_KIND;
  }
  std::memset(out, 0, sizeof(*out));
  out->kind = eng->kind;
  out->version_hash = eon_relax_version_hash();
  out->surface_epoch = eng->epoch;

  auto pot = std::make_shared<eonc::HostSurfacePotential>(surface, user,
                                                          eng->epoch);
  try {
    if (eng->kind == EON_RELAX_KIND_NEB) {
      if (band->n_images < 3) {
        return EON_RELAX_BAND_TOO_SHORT;
      }
      const long want = eng->params.neb_options.image_count + 2;
      if (band->n_images != want) {
        return EON_RELAX_BAND_SIZE;
      }
      std::vector<Matter> path;
      path.reserve(static_cast<size_t>(band->n_images));
      for (long i = 0; i < band->n_images; ++i) {
        Matter img(pot, eng->params);
        fill_matter(img, band, i, eng->epoch);
        path.push_back(std::move(img));
      }
      auto neb = std::make_unique<NudgedElasticBand>(std::move(path),
                                                     eng->params, pot);
      pot->bindPath(neb->path);
      const auto st = neb->compute();
      if (st == NudgedElasticBand::NEBStatus::INIT ||
          st == NudgedElasticBand::NEBStatus::RUNNING) {
        return EON_RELAX_UNKNOWN_STATUS;
      }
      out->status = neb_status(st);
      out->climbing_image = neb->climbingImage;
      out->max_force = neb->convergenceForce();
      out->surface_epoch = pot->surfaceEpoch();
      for (long i = 0; i < band->n_images; ++i) {
        write_image(band->positions + i * 3 * band->n_atoms, *neb->path[i],
                    band->n_atoms);
      }
      return EON_RELAX_OK;
    }

    if (band->n_images != 1) {
      return EON_RELAX_SADDLE_NIMAGES;
    }
    if (!band->mode) {
      return EON_RELAX_NULL_MODE;
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
    out->status = sst;
    return EON_RELAX_OK;
  } catch (const std::exception &) {
    out->status = 0;
    return EON_RELAX_SURFACE_FATAL;
  }
}

void eon_relax_destroy(EonRelaxEngine *eng) {
  delete eng;
}

const char *eon_relax_status_name(eon_relax_kind_t kind, int status) {
  if (kind == EON_RELAX_KIND_NEB) {
    switch (status) {
    case EON_RELAX_NEB_GOOD:
      return "GOOD";
    case EON_RELAX_NEB_INIT:
      return "INIT";
    case EON_RELAX_NEB_BAD_MAX_ITERATIONS:
      return "BAD_MAX_ITERATIONS";
    case EON_RELAX_NEB_RUNNING:
      return "RUNNING";
    case EON_RELAX_NEB_MAX_UNCERTAINTY:
      return "MAX_UNCERTAINTY";
    default:
      return nullptr;
    }
  }
  if (kind == EON_RELAX_KIND_SADDLE) {
    const auto msg = MinModeSaddleSearch::statusMessage(status);
    if (msg == "Unknown status") {
      return nullptr;
    }
    return msg.data();
  }
  return nullptr;
}

} // extern "C"
