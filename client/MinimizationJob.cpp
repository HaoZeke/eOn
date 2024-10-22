#include "MinimizationJob.h"
#include "BaseStructures.h"
#include "HelperFunctions.h"
#include "Matter.h"
#include "Optimizer.h"

std::vector<std::string> MinimizationJob::run(void) {
  std::string posInFilename("pos.con");
  std::string posOutFilename("min.con");

  if (params->main.checkpoint) {
    FILE *pos_file;
    pos_file = fopen("pos_cp.con", "r");
    if (pos_file != NULL) {
      posInFilename = "pos_cp.con";
      SPDLOG_LOGGER_DEBUG(log, "[Minimization] Resuming from checkpoint");
    } else {
      SPDLOG_LOGGER_DEBUG(log, "[Minimization] No checkpoint files found");
    }
  }

  std::vector<std::string> returnFiles;
  returnFiles.push_back(posOutFilename);

  auto pos = std::make_shared<Matter>(pot, params);
  pos->con2matter(posInFilename);

  SPDLOG_LOGGER_DEBUG(log, "\nBeginning minimization of {}", posInFilename);

  bool converged;
  try {
    converged = pos->relax(false, params->debug.writeMovies,
                           params->main.checkpoint, "minimization", "pos");
    if (converged) {
      status = RunStatus::GOOD;
      SPDLOG_LOGGER_DEBUG(log, "Minimization converged within tolerence");
    } else {
      status = RunStatus::FAIL_MAX_ITERATIONS;
      SPDLOG_LOGGER_DEBUG(log, "Minimization did not converge to tolerence!"
                               "Maybe try to increase max_iterations?");
    }
  } catch (int e) {
    if (e == 100) {
      status = RunStatus::FAIL_POTENTIAL_FAILED;
    } else {
      throw e;
    }
  }

  SPDLOG_LOGGER_DEBUG(log, "Saving result to {}", posOutFilename);
  pos->matter2con(posOutFilename);
  if (status != RunStatus::FAIL_POTENTIAL_FAILED) {
    SPDLOG_LOGGER_DEBUG(log, "Final Energy: {}", pos->getPotentialEnergy());
  }

  FILE *fileResults;

  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);
  fileResults = fopen(resultsFilename.c_str(), "wb");

  fprintf(fileResults, "%s termination_reason\n",
          (std::string{magic_enum::enum_name<RunStatus>(status)}).c_str());
  fprintf(fileResults, "minimization job_type\n");
  fprintf(fileResults, "%s potential_type\n",
          std::string{magic_enum::enum_name<PotType>(params->pot.potential)}
              .c_str());
  // fprintf(fileResults, "%d total_force_calls\n", Potential::fcallsTotal);
  if (status != RunStatus::FAIL_POTENTIAL_FAILED) {
    fprintf(fileResults, "%f potential_energy\n", pos->getPotentialEnergy());
  }
  fclose(fileResults);

  return returnFiles;
}
