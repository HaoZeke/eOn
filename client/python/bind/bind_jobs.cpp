#include "Job.h"
#include "Parameters.h"

#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace eonc::pybind {
namespace nb = nanobind;

void bind_jobs(nb::module_ &m) {
  nb::class_<eonc::Job>(m, "Job", "Abstract eOn client job")
      .def("get_type", &eonc::Job::getType)
      .def(
          "run",
          [](eonc::Job &self) { return self.run(); },
          "Run the job in the *current working directory* (config.ini, "
          "pos.con, …). Returns list of output filenames produced.");

  m.def(
      "make_job",
      [](eonc::Parameters &params) {
        auto job =
            eonc::helpers::makeJob(std::make_unique<eonc::Parameters>(params));
        if (!job)
          throw std::runtime_error("make_job: unknown or unsupported job type");
        return std::shared_ptr<eonc::Job>(std::move(job));
      },
      nb::arg("parameters"),
      "Construct a Job from Parameters.main.job (copies Parameters).");

  m.def(
      "run_job_in_directory",
      [](const std::string &workdir, eonc::Parameters params) {
        namespace fs = std::filesystem;
        const fs::path dir(workdir);
        if (!fs::is_directory(dir))
          throw std::invalid_argument("not a directory: " + workdir);

        const fs::path prev = fs::current_path();
        fs::current_path(dir);
        try {
          const fs::path ini = dir / "config.ini";
          if (fs::is_regular_file(ini)) {
            if (params.load(ini.string()) != 0)
              throw std::runtime_error("failed to load " + ini.string());
          }
          auto job =
              eonc::helpers::makeJob(std::make_unique<eonc::Parameters>(params));
          if (!job)
            throw std::runtime_error("unknown job type for makeJob");
          auto files = job->run();
          fs::current_path(prev);
          return files;
        } catch (...) {
          fs::current_path(prev);
          throw;
        }
      },
      nb::arg("workdir"), nb::arg("parameters"),
      "chdir to workdir, optionally reload config.ini, makeJob, run, restore "
      "cwd. Returns output filenames.");
}

} // namespace eonc::pybind
