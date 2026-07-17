#pragma once
/** Shared Matter clone / inplace helpers for pyeonclient bindings.
 *
 * Contract: algorithms never mutate the caller's Matter unless
 * inplace=True. Default is always non-mutating (copy-then-run).
 */
#include "Matter.h"
#include <memory>

namespace eonc::pybind {

/// Working Matter: same object if inplace, otherwise deep copy.
inline std::shared_ptr<eonc::Matter>
matter_work(std::shared_ptr<eonc::Matter> matter, bool inplace) {
  if (inplace)
    return matter;
  return std::make_shared<eonc::Matter>(*matter);
}

/// Fold algorithm result into the caller's Matter when inplace=True.
/// When the job returns a distinct Matter (e.g. ParallelReplica trajectory),
/// copy it into the original shared_ptr so the Python object is updated.
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

} // namespace eonc::pybind
