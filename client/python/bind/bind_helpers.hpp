#pragma once
/** Shared Matter clone / inplace helpers for pyeonclient bindings. */
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

} // namespace eonc::pybind
