#include "Parameters.h"
#include "Job.h"
#include "ProcessSearchJob.h"
#include "SaddleSearchJob.h"
#include "MinimizationJob.h"
#include "PointJob.h"
#include "HessianJob.h"
#include "ParallelReplicaJob.h"
#include "SafeHyperJob.h"
#include "TADJob.h"
#include "ReplicaExchangeJob.h"
#include "BasinHoppingJob.h"
#include "FiniteDifferenceJob.h"
#include "NudgedElasticBandJob.h"
#include "DynamicsJob.h"
#include "PrefactorJob.h"
#include "TestJob.h"
#include "GlobalOptimizationJob.h"
#include "StructureComparisonJob.h"
#include "MonteCarloJob.h"

namespace helper_functions {
std::unique_ptr<Job> makeJob(std::unique_ptr<Parameters> params) {
    switch (params->job) {
      case JobType::ProcessSearch: {
          return (std::make_unique<ProcessSearchJob>(std::move(params)));
          break;
      }
      case JobType::SaddleSearch: {
          return (std::make_unique<SaddleSearchJob>(std::move(params)));
          break;
      }
      case JobType::Minimization: {
          return (std::make_unique<MinimizationJob>(std::move(params)));
          break;
      }
      case JobType::Point: {
          return (std::make_unique<PointJob>(std::move(params)));
          break;
      }
      case JobType::ParallelReplica: {
          return (std::make_unique<ParallelReplicaJob>(std::move(params)));
          break;
      }
      case JobType::SafeHyperdynamics: {
          return (std::make_unique<SafeHyperJob>(std::move(params)));
          break;
      }
      case JobType::TAD: {
          return (std::make_unique<TADJob>(std::move(params)));
          break;
      }
      case JobType::ReplicaExchange: {
          return (std::make_unique<ReplicaExchangeJob>(std::move(params)));
          break;
      }
      case JobType::BasinHopping: {
          return (std::make_unique<BasinHoppingJob>(std::move(params)));
          break;
      }
      case JobType::Hessian: {
          return (std::make_unique<HessianJob>(std::move(params)));
          break;
      }
      case JobType::FiniteDifference: {
          return (std::make_unique<FiniteDifferenceJob>(std::move(params)));
          break;
      }
      case JobType::NEB: {
          return (std::make_unique<NudgedElasticBandJob>(std::move(params)));
          break;
      }
      case JobType::Dynamics: {
          return (std::make_unique<DynamicsJob>(std::move(params)));
          break;
      }
      case JobType::Prefactor: {
          return (std::make_unique<PrefactorJob>(std::move(params)));
          break;
      }
      case JobType::GlobalOptimization: {
          return (std::make_unique<GlobalOptimizationJob>(std::move(params)));
          break;
      }
      case JobType::StructureComparison: {
          return (std::make_unique<StructureComparisonJob>(std::move(params)));
          break;
      }
      case JobType::MonteCarlo: {
          return (std::make_unique<MonteCarloJob>(std::move(params)));
          break;
      }
      default:
          throw std::runtime_error("No known job could be constructed");
        break;
    }
  }

} // namespace helper_functions
