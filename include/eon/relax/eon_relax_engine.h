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
 * - Closed enums. Unknown kind is EON_RELAX_UNKNOWN_KIND. Other
 *   unknown enumerants are EON_RELAX_INVALID_PARAMETER. Never a
 *   silent fallback.
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

#define EON_RELAX_ABI_MAJOR 1u
#define EON_RELAX_ABI_MINOR 1u
#define EON_RELAX_ABI_LAYOUT 1u
#define EON_RELAX_ABI_VERSION                                                  \
  ((EON_RELAX_ABI_MAJOR << 24) | (EON_RELAX_ABI_MINOR << 16) |                 \
   (EON_RELAX_ABI_LAYOUT))
#define EON_RELAX_ABI_UNPACK_MAJOR(v) ((uint32_t)(((v) >> 24) & 0xffu))
#define EON_RELAX_ABI_UNPACK_MINOR(v) ((uint32_t)(((v) >> 16) & 0xffu))
#define EON_RELAX_ABI_UNPACK_LAYOUT(v) ((uint32_t)((v) & 0xffffu))

typedef struct {
  uint32_t major;
  uint32_t minor;
  uint32_t layout_revision;
} eon_relax_abi_stamp_t;

typedef enum {
  EON_RELAX_KIND_INVALID = 0,
  EON_RELAX_KIND_NEB = 1,
  EON_RELAX_KIND_SADDLE = 2
} eon_relax_kind_t;

#define EON_RELAX_KIND_IS_KNOWN(k)                                             \
  ((k) == EON_RELAX_KIND_NEB || (k) == EON_RELAX_KIND_SADDLE)

/** Call rc: 0 success (read outcome), positive recoverable, negative fatal. */
typedef enum {
  EON_RELAX_OK = 0,
  EON_RELAX_NULL_ENGINE = -1,
  EON_RELAX_NULL_BAND = -2,
  EON_RELAX_NULL_SURFACE = -3,
  EON_RELAX_NULL_OUTCOME = -4,
  EON_RELAX_NULL_POSITIONS = -5,
  EON_RELAX_NULL_ATOMIC_NRS = -6,
  EON_RELAX_NULL_BOXES = -7,
  EON_RELAX_NULL_MODE = -8,
  EON_RELAX_BAND_TOO_SHORT = -9,
  EON_RELAX_NATOMS = -10,
  EON_RELAX_UNKNOWN_KIND = -11,
  EON_RELAX_UNKNOWN_STATUS = -12,
  EON_RELAX_CAPNP_ROOT = -13,
  EON_RELAX_SURFACE_FATAL = -15,
  EON_RELAX_ALLOC = -19,
  EON_RELAX_SADDLE_NIMAGES = -17,
  EON_RELAX_BAND_SIZE = -18,
  EON_RELAX_UNAVAILABLE = -21,
  EON_RELAX_INVALID_PARAMETER = -22
} eon_relax_rc_t;

/** NudgedElasticBand::NEBStatus integers. Written to outcome.status for NEB. */
typedef enum {
  EON_RELAX_NEB_GOOD = 0,
  EON_RELAX_NEB_INIT = 1,
  EON_RELAX_NEB_BAD_MAX_ITERATIONS = 2,
  EON_RELAX_NEB_RUNNING = 3,
  EON_RELAX_NEB_MAX_UNCERTAINTY = 4
} eon_relax_neb_status_t;

#define EON_RELAX_IMAGE_NONE (-1L)

/**
 * Caller-owned band. positions is 3*n_atoms*n_images (image stride
 * 3*n_atoms) and is written on a successful run. atomic_nrs is n_atoms
 * (shared composition), boxes is 9*n_images (row-major). is_fixed may
 * be NULL (all free) or n_atoms (1 = fixed atom). mode may be NULL;
 * required for kind=saddle (3*n_atoms). image_ids may be NULL.
 * positions layout matches Matter AtomMatrix (row-major xyz per atom).
 * boxes are 9 doubles per image, same order as Matter::getCell().data().
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
  eon_relax_kind_t kind;
  int status;
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
EON_RELAX_API int eon_relax_available(void);
EON_RELAX_API int eon_relax_abi_stamp(eon_relax_abi_stamp_t *out);
/** Immortal dest identity string. Never NULL. */
EON_RELAX_API const char *eon_relax_version_hash_str(void);
/** FNV-1a64 of eon_relax_version_hash_str(). Copied to outcome.version_hash. */
EON_RELAX_API uint64_t eon_relax_version_hash(void);

/**
 * Create from a Cap'n Proto flat-array RelaxEngineParams message.
 * config may be NULL (documented defaults, kind=NEB). Returns NULL on
 * failure with a message in errbuf.
 */
EON_RELAX_API EonRelaxEngine *eon_relax_create(const void *config,
                                               size_t config_len, char *errbuf,
                                               size_t errlen);

EON_RELAX_API int eon_relax_set_surface_epoch(EonRelaxEngine *eng,
                                              uint64_t epoch);
/** Last named run/create rc on this handle. 0 if none. */
EON_RELAX_API int eon_relax_last_error(const EonRelaxEngine *eng);

/**
 * Run the bound engine. Writes updated coordinates into
 * band->positions (caller-owned). surface may not be NULL. Kind is
 * bound at create (NULL config => NEB). Returns eon_relax_rc_t;
 * on 0, out->status is the dest NEB or saddle table.
 */
EON_RELAX_API int eon_relax_run(EonRelaxEngine *eng, eon_relax_band_t *band,
                                eon_relax_surface_fn surface, void *user,
                                eon_relax_outcome_t *out);

EON_RELAX_API void eon_relax_destroy(EonRelaxEngine *eng);

/**
 * One band optimizer step (kind=NEB only; minor >= 1). The first call
 * on a fresh or reset engine initializes the band state from the
 * caller's positions and captures the baseline force; every call
 * syncs band->positions in (the host may have moved images between
 * steps), assembles NEB forces on the host surface, takes ONE
 * optimizer step, and writes the stepped positions back. Climbing
 * image activation follows the engine's own trigger rule; saddle
 * POLICY (MMF bursts, resampling, acquisition) stays with the host.
 * outcome.status is RUNNING until convergenceForce() reaches the
 * configured tolerance, then GOOD. surface and user must not change
 * across steps without an eon_relax_reset.
 */
EON_RELAX_API int eon_relax_step(EonRelaxEngine *eng, eon_relax_band_t *band,
                                 eon_relax_surface_fn surface, void *user,
                                 eon_relax_outcome_t *out);

/**
 * Drop the stepper's band state and optimizer history; the next
 * eon_relax_step reinitializes from the caller's band. This is the
 * host's model-update hook: quasi-Newton memory taken on one surface
 * epoch must not survive onto the next.
 */
EON_RELAX_API int eon_relax_reset(EonRelaxEngine *eng);

/** NULL if kind or status is unknown. */
EON_RELAX_API const char *eon_relax_status_name(eon_relax_kind_t kind,
                                                int status);

#ifdef __cplusplus
}
#endif
