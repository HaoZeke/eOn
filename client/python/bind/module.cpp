/*
** pyeonclient._core — complete nanobind surface for the eOn C++ client.
*/
#include <nanobind/nanobind.h>

namespace eonc::pybind {
namespace nb = nanobind;

void bind_enums(nb::module_ &m);
void bind_parameters(nb::module_ &m);
void bind_potential(nb::module_ &m);
void bind_matter(nb::module_ &m);
void bind_jobs(nb::module_ &m);
void bind_neb(nb::module_ &m);

} // namespace eonc::pybind

NB_MODULE(_core, m) {
  m.doc() =
      "eOn client core via nanobind: Matter, Parameters, Potential, Jobs. "
      "Stable ABI (abi3) or free-threaded; ConFileIO/readcon for .con I/O.";
  m.attr("__version__") = "0.2.0";

  eonc::pybind::bind_enums(m);
  eonc::pybind::bind_parameters(m);
  eonc::pybind::bind_potential(m);
  eonc::pybind::bind_matter(m);
  eonc::pybind::bind_jobs(m);
  eonc::pybind::bind_neb(m);
}
