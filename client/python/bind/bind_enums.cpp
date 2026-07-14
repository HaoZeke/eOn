#include "BaseStructures.h"
#include "ConFileIO.h"
#include "Matter.h"

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <magic_enum/magic_enum.hpp>
#include <string>

namespace eonc::pybind {
namespace nb = nanobind;

void bind_enums(nb::module_ &m) {
  nb::enum_<eonc::io::IoStatus>(m, "IoStatus",
                               "I/O result from ConFileIO / Matter file ops")
      .value("Ok", eonc::io::IoStatus::Ok)
      .value("ReadError", eonc::io::IoStatus::ReadError)
      .value("WriteError", eonc::io::IoStatus::WriteError)
      .value("AppendError", eonc::io::IoStatus::AppendError)
      .value("OpenError", eonc::io::IoStatus::OpenError)
      .value("InvalidArgument", eonc::io::IoStatus::InvalidArgument)
      .export_values();

  m.def(
      "io_ok",
      [](eonc::io::IoStatus s) { return eonc::io::io_ok(s); },
      nb::arg("status"), "True if status is IoStatus.Ok");

  m.def(
      "io_status_name",
      [](eonc::io::IoStatus s) {
        return std::string(eonc::io::io_status_name(s));
      },
      nb::arg("status"));

  nb::enum_<eonc::PbcConvention>(m, "PbcConvention",
                                 "Position wrap convention (issue #176)")
      .value("Legacy", eonc::PbcConvention::Legacy)
      .value("MinimumImage", eonc::PbcConvention::MinimumImage)
      .export_values();

  nb::enum_<eonc::PotType>(m, "PotType", "Potential type (matches config names)")
      .value("UNKNOWN", eonc::PotType::UNKNOWN)
      .value("EMT", eonc::PotType::EMT)
      .value("EXT_POT", eonc::PotType::EXT_POT)
      .value("LJ", eonc::PotType::LJ)
      .value("LJCLUSTER", eonc::PotType::LJCLUSTER)
      .value("MORSE_PT", eonc::PotType::MORSE_PT)
      .value("CUH2", eonc::PotType::CUH2)
      .value("EAM_AL", eonc::PotType::EAM_AL)
      .value("EDIP", eonc::PotType::EDIP)
      .value("FEHE", eonc::PotType::FEHE)
      .value("LENOSKY_SI", eonc::PotType::LENOSKY_SI)
      .value("SW_SI", eonc::PotType::SW_SI)
      .value("TERSOFF_SI", eonc::PotType::TERSOFF_SI)
      .value("LAMMPS", eonc::PotType::LAMMPS)
      .value("XTB", eonc::PotType::XTB)
      .value("ASE_ORCA", eonc::PotType::ASE_ORCA)
      .value("ASE_POT", eonc::PotType::ASE_POT)
      .value("ASE_NWCHEM", eonc::PotType::ASE_NWCHEM)
      .value("METATOMIC", eonc::PotType::METATOMIC)
      .value("ZBL", eonc::PotType::ZBL)
      .value("SocketNWChem", eonc::PotType::SocketNWChem)
      .value("RGPOT", eonc::PotType::RGPOT)
      .value("CatLearn", eonc::PotType::CatLearn)
      .value("QSC", eonc::PotType::QSC)
      .export_values();

  m.def(
      "pot_type_from_name",
      [](const std::string &name) {
        auto v = magic_enum::enum_cast<eonc::PotType>(
            name, magic_enum::case_insensitive);
        if (!v) {
          throw std::invalid_argument("unknown PotType: " + name);
        }
        return *v;
      },
      nb::arg("name"), "Parse PotType from config-style name (case-insensitive)");

  m.def(
      "pot_type_name",
      [](eonc::PotType t) {
        return std::string(magic_enum::enum_name(t));
      },
      nb::arg("pot_type"));

  nb::enum_<eonc::JobType>(m, "JobType")
      .value("Unknown", eonc::JobType::Unknown)
      .value("Process_Search", eonc::JobType::Process_Search)
      .value("Saddle_Search", eonc::JobType::Saddle_Search)
      .value("Minimization", eonc::JobType::Minimization)
      .value("Point", eonc::JobType::Point)
      .value("Parallel_Replica", eonc::JobType::Parallel_Replica)
      .value("Basin_Hopping", eonc::JobType::Basin_Hopping)
      .value("Hessian", eonc::JobType::Hessian)
      .value("Nudged_Elastic_Band", eonc::JobType::Nudged_Elastic_Band)
      .value("Dynamics", eonc::JobType::Dynamics)
      .value("Prefactor", eonc::JobType::Prefactor)
      .value("Monte_Carlo", eonc::JobType::Monte_Carlo)
      .export_values();

  nb::enum_<eonc::RunStatus>(m, "RunStatus")
      .value("GOOD", eonc::RunStatus::GOOD)
      .value("FAIL_MAX_ITERATIONS", eonc::RunStatus::FAIL_MAX_ITERATIONS)
      .value("FAIL_POTENTIAL_FAILED", eonc::RunStatus::FAIL_POTENTIAL_FAILED)
      .export_values();
}

} // namespace eonc::pybind
