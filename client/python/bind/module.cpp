/*
** pyeonclient._core — nanobind bindings for eOn C++ client internals.
**
** Bottom-up surface: enums, Parameters, Potential, Matter.
** Jobs / communicator adapters layer on top of Matter (later PRs).
*/
#include <nanobind/nanobind.h>

namespace eonc::pybind {
namespace nb = nanobind;

void bind_enums(nb::module_ &m);
void bind_parameters(nb::module_ &m);
void bind_potential(nb::module_ &m);
void bind_matter(nb::module_ &m);

} // namespace eonc::pybind

NB_MODULE(_core, m) {
  m.doc() =
      "eOn client core (Matter, Parameters, Potential) via nanobind. "
      "I/O uses ConFileIO / readcon-core — the same path as eonclient.";
  m.attr("__version__") = "0.1.0";

  eonc::pybind::bind_enums(m);
  eonc::pybind::bind_parameters(m);
  eonc::pybind::bind_potential(m);
  eonc::pybind::bind_matter(m);
}
