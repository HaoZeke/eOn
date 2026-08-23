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
/**
 * @file eon_relax_engine.h
 * @brief Tier-2 C ABI: eOn as a dlopen relaxation engine.
 *
 * Reverse of include/eon/potentials/Rgpot/engine_c_abi.h: the host
 * supplies the surface (energy/forces/variance for a band), eOn owns
 * NudgedElasticBand / MinModeSaddleSearch. Configuration is a Cap'n
 * Proto flat-array whose root is RelaxEngineParams in
 * schema/eon_params.capnp (NebParams / SaddleParams arms). NULL config
 * selects Parameters NSDMI defaults and kind=NEB.
 *
 * Units: positions Angstrom, energy eV, forces eV/Angstrom.
 *
 * Conventions:
 *
 * - Return values are three-valued: 0 success, positive recoverable,
 *   negative fatal (the instance is unusable).
 * - Closed enums. Unknown enumerants are fail-closed
 *   (EON_RELAX_INVALID_PARAMETER), never a silent fallback.
 * - Version stamp: major = layout break (do not read other fields),
 *   minor = new enum values / additive symbols.
 * - Capability discovery is by symbol presence (dlsym).
 */
#pragma once

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#ifdef EON_RELAX_ENGINE_BUILD
#define EON_RELAX_API __declspec(dllexport)
#else
#define EON_RELAX_API __declspec(dllimport)
#endif
#else
#define EON_RELAX_API __attribute__((visibility("default")))
#endif

typedef struct EonRelaxEngine EonRelaxEngine;

#define EON_RELAX_ABI_MAJOR 1
#define EON_RELAX_ABI_MINOR 0
#define EON_RELAX_ABI_LAYOUT 1
#define EON_RELAX_ABI_VERSION                                                  \
  ((EON_RELAX_ABI_MAJOR << 16) | (EON_RELAX_ABI_MINOR))

typedef struct {
  uint32_t major;
  uint32_t minor;
  uint32_t layout_revision;
} eon_relax_abi_stamp_t;

typedef enum {
  EON_RELAX_KIND_NEB = 0,
  EON_RELAX_KIND_SADDLE = 1
} eon_relax_kind_t;

typedef enum {
  EON_RELAX_OK = 0,
  EON_RELAX_CONVERGED = 1,
  EON_RELAX_MAX_ITERATIONS = 2,
  EON_RELAX_MAX_UNCERTAINTY = 3,
  EON_RELAX_INIT = 4,
  EON_RELAX_INVALID_PARAMETER = -1,
  EON_RELAX_UNAVAILABLE = -2,
  EON_RELAX_SURFACE_FAILED = -3,
  EON_RELAX_UNKNOWN_KIND = -4
} eon_relax_status_t;

/**
 * Caller-owned band. positions is 3*n_atoms*n_images (image stride
 * 3*n_atoms) and is written on a successful run. atomic_nrs is n_atoms
 * (shared composition), boxes is 9*n_images (row-major). is_fixed may
 * be NULL (all free) or n_atoms (1 = fixed atom). mode may be NULL;
 * required for kind=saddle (3*n_atoms). image_ids may be NULL
 * (then 0..n_images-1).
 */
typedef struct {
  long n_images;
  long n_atoms;
  double *positions;
  const int *atomic_nrs;
  const double *boxes;
  const int *is_fixed;
  const double *mode;
  const long *image_ids;
} eon_relax_band_t;

typedef struct {
  eon_relax_status_t status;
  long iterations;
  long climbing_image;
  double max_force;
  uint64_t version_hash;
  uint64_t surface_epoch;
} eon_relax_outcome_t;

/**
 * Host surface. n_images systems, one composition. positions/forces
 * 3*n_atoms*n_images, boxes 9*n_images, energies n_images, variances
 * n_images (may be NULL). image_ids n_images (may be NULL). epoch_out
 * may be NULL; if written, Matter caches key on the new generation.
 * Return 0 on success, positive recoverable, negative fatal.
 */
typedef int (*eon_relax_surface_fn)(
    void *user, long n_images, long n_atoms, const double *positions,
    const int *atomic_nrs, const double *boxes, const long *image_ids,
    double *energies, double *forces, double *variances, uint64_t *epoch_out);

EON_RELAX_API int eon_relax_abi_version(void);
EON_RELAX_API int eon_relax_abi_stamp(eon_relax_abi_stamp_t *out);
EON_RELAX_API uint64_t eon_relax_version_hash(void);

/**
 * Create from a Cap'n Proto flat-array RelaxEngineParams message.
 * config may be NULL (documented defaults, kind=NEB). Returns NULL on
 * failure with a message in errbuf.
 */
EON_RELAX_API EonRelaxEngine *eon_relax_create(const void *config,
                                               size_t config_len, char *errbuf,
                                               size_t errlen);

/**
 * Select NEB or saddle. Unknown kind is fail-closed. Default after
 * NULL create is NEB.
 */
EON_RELAX_API int eon_relax_set_kind(EonRelaxEngine *eng,
                                     eon_relax_kind_t kind);

EON_RELAX_API int eon_relax_set_surface_epoch(EonRelaxEngine *eng,
                                              uint64_t epoch);

/**
 * Run the bound engine. Writes updated coordinates into
 * band->positions (caller-owned). surface may not be NULL.
 */
EON_RELAX_API int eon_relax_run(EonRelaxEngine *eng,
                                const eon_relax_band_t *band,
                                eon_relax_surface_fn surface, void *user,
                                eon_relax_outcome_t *out);

EON_RELAX_API void eon_relax_destroy(EonRelaxEngine *eng);

#ifdef __cplusplus
}
#endif
