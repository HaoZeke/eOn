#pragma once
/** Shared Matter helpers for pyeonclient bindings only (not on Matter itself).
 *
 * Contract: algorithms never mutate the caller's Matter unless
 * inplace=True. Default is always non-mutating (copy-then-run).
 *
 * Buffer access uses public Matter getters/setters; ASE/Python cache state
 * lives in the binding layer (AseCalcPotential), not on Matter.
 */
#include "Eigen.h"
#include "Matter.h"

#include <memory>
#include <stdexcept>

namespace eonc::pybind {

/// Working Matter: same object if inplace, otherwise deep copy.
inline std::shared_ptr<eonc::Matter>
matter_work(std::shared_ptr<eonc::Matter> matter, bool inplace) {
  if (inplace)
    return matter;
  return std::make_shared<eonc::Matter>(*matter);
}

/// Fold algorithm result into the caller's Matter when inplace=True.
inline std::shared_ptr<eonc::Matter>
matter_result(std::shared_ptr<eonc::Matter> original,
              std::shared_ptr<eonc::Matter> result, bool inplace) {
  if (!result)
    return result;
  if (!inplace)
    return result;
  if (result.get() != original.get())
    *original = *result;
  return original;
}

// Binding-local zero-copy views / bulk sets via public Matter API.

inline double *matter_positions_ptr(eonc::Matter &m) {
  return const_cast<double *>(m.getPositions().data());
}
inline const double *matter_positions_ptr(const eonc::Matter &m) {
  return m.getPositions().data();
}
inline double *matter_forces_ptr(eonc::Matter &m) {
  return const_cast<double *>(m.getForces().data());
}
inline double *matter_forces_raw_ptr(eonc::Matter &m) {
  return const_cast<double *>(m.getForcesRaw().data());
}

inline void matter_set_positions_buf(eonc::Matter &m, const double *xyz,
                                     long n) {
  if (n != m.numberOfAtoms()) {
    throw std::invalid_argument("positions n != n_atoms");
  }
  m.setPositions(AtomMatrix::Map(const_cast<double *>(xyz), n, 3));
}

inline void matter_set_cell_buf(eonc::Matter &m, const double *c33) {
  m.setCell(Matrix3d::Map(const_cast<double *>(c33)));
}

inline void matter_set_velocities_buf(eonc::Matter &m, const double *v,
                                      long n) {
  if (n != m.numberOfAtoms()) {
    throw std::invalid_argument("velocities n != n_atoms");
  }
  m.setVelocities(AtomMatrix::Map(const_cast<double *>(v), n, 3));
}

inline void matter_set_masses_buf(eonc::Matter &m, const double *mass, long n) {
  if (n != m.numberOfAtoms()) {
    throw std::invalid_argument("masses n != n_atoms");
  }
  m.setMasses(
      Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1>>(mass, n));
}

inline void matter_set_z_buf(eonc::Matter &m, const int *z, long n) {
  if (n != m.numberOfAtoms()) {
    throw std::invalid_argument("atomic_numbers n != n_atoms");
  }
  m.setAtomicNrs(VectorXi::Map(const_cast<int *>(z), n));
}

inline void matter_set_fixed_buf(eonc::Matter &m, const int *fx, long n) {
  if (n != m.numberOfAtoms()) {
    throw std::invalid_argument("fixed n != n_atoms");
  }
  for (long i = 0; i < n; ++i) {
    m.setFixed(i, fx[i] ? 1 : 0);
  }
}

} // namespace eonc::pybind
