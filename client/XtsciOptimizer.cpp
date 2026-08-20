/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
*/
#include "eon/XtsciOptimizer.h"

#include <algorithm>
#include <stdexcept>

#include <xts.h>

#if XTS_ABI_LAYOUT_REVISION < 2
#error "xtsci-optimize ABI layout 2 or newer is required for max_move"
#endif

namespace {

struct XtsObjectiveContext {
  eonc::ObjectiveFunction *objective;
};

Eigen::Map<const Eigen::VectorXd>
map_input(const DLManagedTensorVersioned *tensor) {
  const auto &dl = tensor->dl_tensor;
  if (dl.ndim != 1 || dl.dtype.code != kDLFloat || dl.dtype.bits != 64 ||
      dl.dtype.lanes != 1 || dl.device.device_type != kDLCPU ||
      dl.shape == nullptr || dl.data == nullptr) {
    throw std::runtime_error("xtsci objective requires a CPU f64 vector");
  }
  return {static_cast<const double *>(dl.data) + dl.byte_offset / sizeof(double),
          static_cast<Eigen::Index>(dl.shape[0])};
}

Eigen::Map<Eigen::VectorXd>
map_output(DLManagedTensorVersioned *tensor) {
  const auto &dl = tensor->dl_tensor;
  if (dl.ndim != 1 || dl.dtype.code != kDLFloat || dl.dtype.bits != 64 ||
      dl.dtype.lanes != 1 || dl.device.device_type != kDLCPU ||
      dl.shape == nullptr || dl.data == nullptr) {
    throw std::runtime_error("xtsci gradient requires a CPU f64 vector");
  }
  return {static_cast<double *>(dl.data) + dl.byte_offset / sizeof(double),
          static_cast<Eigen::Index>(dl.shape[0])};
}

xts_status_t evaluate(void *user, const DLManagedTensorVersioned *x,
                      double *value_out) {
  try {
    auto *context = static_cast<XtsObjectiveContext *>(user);
    const auto positions = map_input(x);
    context->objective->setPositions(positions);
    *value_out = context->objective->getEnergy();
    return XTS_SUCCESS;
  } catch (...) {
    return XTS_INTERNAL_ERROR;
  }
}

xts_status_t gradient(void *user, const DLManagedTensorVersioned *x,
                      DLManagedTensorVersioned *gradient_out) {
  try {
    auto *context = static_cast<XtsObjectiveContext *>(user);
    const auto positions = map_input(x);
    context->objective->setPositions(positions);
    const auto gradient_value = context->objective->getGradient();
    auto output = map_output(gradient_out);
    if (output.size() != gradient_value.size()) {
      return XTS_INVALID_PARAMETER;
    }
    output = gradient_value;
    return XTS_SUCCESS;
  } catch (...) {
    return XTS_INTERNAL_ERROR;
  }
}

} // namespace

namespace eonc {

int XtsciOptimizer::step(double a_maxMove) {
  return run(1, a_maxMove);
}

int XtsciOptimizer::run(size_t a_maxIterations, double a_maxMove) {
  if (m_objf->degreesOfFreedom() <= 0 || a_maxIterations == 0) {
    return m_objf->isConverged() ? 1 : 0;
  }

  auto positions = m_objf->getPositions();
  if (positions.size() != m_objf->degreesOfFreedom()) {
    throw std::runtime_error("xtsci objective position dimension mismatch");
  }
  const auto stamp = xts_abi_stamp();
  if (xts_abi_compatible(&stamp) == 0) {
    throw std::runtime_error("incompatible xtsci-optimize ABI");
  }
  auto *tensor = xts_tensor_borrow_cpu_f64(positions.data(), positions.size());
  if (tensor == nullptr) {
    throw std::runtime_error("could not allocate xtsci objective tensor");
  }
  XtsObjectiveContext context{m_objf.get()};
  xts_control_t control{
      a_maxIterations,
      m_optConfig.opts.converged_force,
      std::max(a_maxMove, 1.0e-12),
      static_cast<size_t>(std::max<long>(0, m_optConfig.opts.lbfgs.memory)),
      std::max(a_maxMove, 0.0),
  };
  xts_report_t report{};
  const auto status = xts_minimize(
      evaluate, gradient, &context, tensor, &control, XTS_LBFGS, &report);
  xts_tensor_free(tensor);
  if (status != XTS_SUCCESS) {
    throw std::runtime_error(xts_last_error());
  }
  m_objf->setPositions(positions);
  return m_objf->isConverged() ? 1 : 0;
}

} // namespace eonc
