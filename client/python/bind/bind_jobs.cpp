#include "HelperFunctions.h"
#include "Job.h"
#include "Parameters.h"
#include "PotRegistry.h"

#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace eonc::pybind {
namespace nb = nanobind;

namespace {

/** Mirror ClientEON post-job side effects that jobs themselves do not write. */
void append_client_postamble(std::vector<std::string> &files,
                             std::chrono::steady_clock::time_point start) {
  // Potential must be destroyed so PotRegistry records force-call tallies.
  PotRegistry::get().write_summary();
  files.emplace_back("_potcalls.json");

  const auto end = std::chrono::steady_clock::now();
  const std::chrono::duration<double> elapsed = end - start;
  double utime = 0, stime = 0, rtime = 0;
  eonc::helpers::getTime(&rtime, &utime, &stime);

  std::ofstream result_file("results.dat", std::ios::app);
  if (!result_file.is_open())
    throw std::runtime_error("failed to append timing to results.dat");
  result_file << "time_seconds " << elapsed.count() << "\n";
#ifndef _WIN32
  result_file << "user_time " << utime << "\n";
  result_file << "system_time " << stime << "\n";
#endif
}

} // namespace

void bind_jobs(nb::module_ &m) {
  nb::class_<eonc::Job>(m, "Job", "Abstract eOn client job")
      .def("get_type", &eonc::Job::getType)
      .def(
          "run",
          [](eonc::Job &self) { return self.run(); },
          "Run the job in the *current working directory* (config.ini, "
          "pos.con, …). Returns list of output filenames produced. "
          "Does not append ClientEON timing / potcall summary — prefer "
          "run_job_in_directory for eonclient parity.");

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
          const auto start = std::chrono::steady_clock::now();
          auto job =
              eonc::helpers::makeJob(std::make_unique<eonc::Parameters>(params));
          if (!job)
            throw std::runtime_error("unknown job type for makeJob");
          auto files = job->run();
          // Destroy Potential before potcall summary (same order as ClientEON).
          job.reset();
          append_client_postamble(files, start);
          fs::current_path(prev);
          return files;
        } catch (...) {
          fs::current_path(prev);
          throw;
        }
      },
      nb::arg("workdir"), nb::arg("parameters"),
      "chdir to workdir, optionally reload config.ini, makeJob, run, then "
      "append ClientEON-equivalent timing and _potcalls.json. Returns output "
      "filenames.");
}

} // namespace eonc::pybind
