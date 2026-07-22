"""AUTO-GENERATED from schema/eon_params.capnp — do not edit.

Regenerate: python tools/params_ssot/codegen.py
"""
from __future__ import annotations

CATALOG = {
  "flat_aliases": [
    {
      "default": 20,
      "flat_key": "lbfgs_memory",
      "ssot_path": "Optimizer.LBFGS.memory",
      "yaml_section": "Optimizer"
    },
    {
      "default": 0.01,
      "flat_key": "lbfgs_inverse_curvature",
      "ssot_path": "Optimizer.LBFGS.inverse_curvature",
      "yaml_section": "Optimizer"
    },
    {
      "default": 0.0,
      "flat_key": "lbfgs_max_inverse_curvature",
      "ssot_path": "Optimizer.LBFGS.max_inverse_curvature",
      "yaml_section": "Optimizer"
    },
    {
      "default": True,
      "flat_key": "lbfgs_auto_scale",
      "ssot_path": "Optimizer.LBFGS.auto_scale",
      "yaml_section": "Optimizer"
    },
    {
      "default": True,
      "flat_key": "lbfgs_angle_reset",
      "ssot_path": "Optimizer.LBFGS.angle_reset",
      "yaml_section": "Optimizer"
    },
    {
      "default": True,
      "flat_key": "lbfgs_distance_reset",
      "ssot_path": "Optimizer.LBFGS.distance_reset",
      "yaml_section": "Optimizer"
    },
    {
      "default": False,
      "flat_key": "cg_no_overshooting",
      "ssot_path": "Optimizer.CG.no_overshooting",
      "yaml_section": "Optimizer"
    },
    {
      "default": False,
      "flat_key": "cg_knock_out_max_move",
      "ssot_path": "Optimizer.CG.knock_out_max_move",
      "yaml_section": "Optimizer"
    },
    {
      "default": False,
      "flat_key": "cg_line_search",
      "ssot_path": "Optimizer.CG.line_search",
      "yaml_section": "Optimizer"
    },
    {
      "default": 0.1,
      "flat_key": "cg_line_converged",
      "ssot_path": "Optimizer.CG.line_converged",
      "yaml_section": "Optimizer"
    },
    {
      "default": 10,
      "flat_key": "cg_max_iter_line_search",
      "ssot_path": "Optimizer.CG.line_search_max_iter",
      "yaml_section": "Optimizer"
    },
    {
      "default": 0,
      "flat_key": "cg_max_iter_before_reset",
      "ssot_path": "Optimizer.CG.max_iter_before_reset",
      "yaml_section": "Optimizer"
    },
    {
      "default": False,
      "flat_key": "qm_steepest_descent",
      "ssot_path": "Optimizer.Quickmin.steepest_descent",
      "yaml_section": "Optimizer"
    },
    {
      "default": 0.1,
      "flat_key": "sd_alpha",
      "ssot_path": "Optimizer.SD.alpha",
      "yaml_section": "Optimizer"
    },
    {
      "default": False,
      "flat_key": "sd_two_point",
      "ssot_path": "Optimizer.SD.two_point",
      "yaml_section": "Optimizer"
    }
  ],
  "schema_version": 1,
  "sections": {
    "Main": {
      "fields": [
        {
          "capnp": "job",
          "default": "process_search",
          "ordinal": 0,
          "snake": "job",
          "type": "Text"
        },
        {
          "capnp": "randomSeed",
          "default": -1,
          "ordinal": 1,
          "snake": "random_seed",
          "type": "Int64"
        },
        {
          "capnp": "temperature",
          "default": 300.0,
          "ordinal": 2,
          "snake": "temperature",
          "type": "Float64"
        },
        {
          "capnp": "quiet",
          "default": False,
          "ordinal": 3,
          "snake": "quiet",
          "type": "Bool"
        },
        {
          "capnp": "writeLog",
          "default": True,
          "ordinal": 4,
          "snake": "write_log",
          "type": "Bool"
        },
        {
          "capnp": "checkpoint",
          "default": False,
          "ordinal": 5,
          "snake": "checkpoint",
          "type": "Bool"
        },
        {
          "capnp": "iniFilename",
          "default": "config.ini",
          "ordinal": 6,
          "snake": "ini_filename",
          "type": "Text"
        },
        {
          "capnp": "conFilename",
          "default": "pos.con",
          "ordinal": 7,
          "snake": "con_filename",
          "type": "Text"
        },
        {
          "capnp": "finiteDifference",
          "default": 0.01,
          "ordinal": 8,
          "snake": "finite_difference",
          "type": "Float64"
        },
        {
          "capnp": "maxForceCalls",
          "default": 0,
          "ordinal": 9,
          "snake": "max_force_calls",
          "type": "Int64"
        },
        {
          "capnp": "removeNetForce",
          "default": True,
          "ordinal": 10,
          "snake": "remove_net_force",
          "type": "Bool"
        },
        {
          "capnp": "writeConForces",
          "default": False,
          "ordinal": 11,
          "snake": "write_con_forces",
          "type": "Bool"
        },
        {
          "capnp": "parallel",
          "default": True,
          "ordinal": 12,
          "snake": "parallel",
          "type": "Bool"
        }
      ],
      "struct": "MainOptions"
    },
    "Optimizer": {
      "fields": [
        {
          "capnp": "optMethod",
          "default": "cg",
          "ordinal": 0,
          "snake": "opt_method",
          "type": "Text"
        },
        {
          "capnp": "convergenceMetric",
          "default": "norm",
          "ordinal": 1,
          "snake": "convergence_metric",
          "type": "Text"
        },
        {
          "capnp": "maxIterations",
          "default": 1000,
          "ordinal": 2,
          "snake": "max_iterations",
          "type": "UInt64"
        },
        {
          "capnp": "maxMove",
          "default": 0.2,
          "ordinal": 3,
          "snake": "max_move",
          "type": "Float64"
        },
        {
          "capnp": "convergedForce",
          "default": 0.01,
          "ordinal": 4,
          "snake": "converged_force",
          "type": "Float64"
        },
        {
          "capnp": "timeStep",
          "default": 1.0,
          "ordinal": 5,
          "snake": "time_step",
          "type": "Float64"
        },
        {
          "capnp": "maxTimeStep",
          "default": 2.5,
          "ordinal": 6,
          "snake": "max_time_step",
          "type": "Float64"
        },
        {
          "capnp": "lbfgs",
          "default": None,
          "nested": True,
          "ordinal": 7,
          "snake": "lbfgs",
          "type": "OptimizerLbfgsOptions"
        },
        {
          "capnp": "cg",
          "default": None,
          "nested": True,
          "ordinal": 8,
          "snake": "cg",
          "type": "OptimizerCgOptions"
        },
        {
          "capnp": "quickmin",
          "default": None,
          "nested": True,
          "ordinal": 9,
          "snake": "quickmin",
          "type": "OptimizerQuickminOptions"
        },
        {
          "capnp": "sd",
          "default": None,
          "nested": True,
          "ordinal": 10,
          "snake": "sd",
          "type": "OptimizerSdOptions"
        }
      ],
      "struct": "OptimizerOptions"
    },
    "Optimizer.CG": {
      "fields": [
        {
          "capnp": "noOvershooting",
          "default": False,
          "ordinal": 0,
          "snake": "no_overshooting",
          "type": "Bool"
        },
        {
          "capnp": "knockOutMaxMove",
          "default": False,
          "ordinal": 1,
          "snake": "knock_out_max_move",
          "type": "Bool"
        },
        {
          "capnp": "lineSearch",
          "default": False,
          "ordinal": 2,
          "snake": "line_search",
          "type": "Bool"
        },
        {
          "capnp": "lineConverged",
          "default": 0.1,
          "ordinal": 3,
          "snake": "line_converged",
          "type": "Float64"
        },
        {
          "capnp": "lineSearchMaxIter",
          "default": 10,
          "ordinal": 4,
          "snake": "line_search_max_iter",
          "type": "Int64"
        },
        {
          "capnp": "maxIterBeforeReset",
          "default": 0,
          "ordinal": 5,
          "snake": "max_iter_before_reset",
          "type": "Int64"
        }
      ],
      "struct": "OptimizerCgOptions"
    },
    "Optimizer.LBFGS": {
      "fields": [
        {
          "capnp": "memory",
          "default": 20,
          "ordinal": 0,
          "snake": "memory",
          "type": "Int64"
        },
        {
          "capnp": "inverseCurvature",
          "default": 0.01,
          "ordinal": 1,
          "snake": "inverse_curvature",
          "type": "Float64"
        },
        {
          "capnp": "maxInverseCurvature",
          "default": 0.0,
          "ordinal": 2,
          "snake": "max_inverse_curvature",
          "type": "Float64"
        },
        {
          "capnp": "autoScale",
          "default": True,
          "ordinal": 3,
          "snake": "auto_scale",
          "type": "Bool"
        },
        {
          "capnp": "angleReset",
          "default": True,
          "ordinal": 4,
          "snake": "angle_reset",
          "type": "Bool"
        },
        {
          "capnp": "distanceReset",
          "default": True,
          "ordinal": 5,
          "snake": "distance_reset",
          "type": "Bool"
        }
      ],
      "struct": "OptimizerLbfgsOptions"
    },
    "Optimizer.Quickmin": {
      "fields": [
        {
          "capnp": "steepestDescent",
          "default": False,
          "ordinal": 0,
          "snake": "steepest_descent",
          "type": "Bool"
        }
      ],
      "struct": "OptimizerQuickminOptions"
    },
    "Optimizer.SD": {
      "fields": [
        {
          "capnp": "alpha",
          "default": 0.1,
          "ordinal": 0,
          "snake": "alpha",
          "type": "Float64"
        },
        {
          "capnp": "twoPoint",
          "default": False,
          "ordinal": 1,
          "snake": "two_point",
          "type": "Bool"
        }
      ],
      "struct": "OptimizerSdOptions"
    },
    "Potential": {
      "fields": [
        {
          "capnp": "potential",
          "default": "lj",
          "ordinal": 0,
          "snake": "potential",
          "type": "Text"
        },
        {
          "capnp": "mpiPollPeriod",
          "default": 0.25,
          "ordinal": 1,
          "snake": "mpi_poll_period",
          "type": "Float64"
        },
        {
          "capnp": "lammpsLogging",
          "default": False,
          "ordinal": 2,
          "snake": "lammps_logging",
          "type": "Bool"
        },
        {
          "capnp": "lammpsThreads",
          "default": 0,
          "ordinal": 3,
          "snake": "lammps_threads",
          "type": "Int32"
        },
        {
          "capnp": "emtRasmussen",
          "default": False,
          "ordinal": 4,
          "snake": "emt_rasmussen",
          "type": "Bool"
        },
        {
          "capnp": "logPotential",
          "default": False,
          "ordinal": 5,
          "snake": "log_potential",
          "type": "Bool"
        },
        {
          "capnp": "extPotPath",
          "default": "./ext_pot",
          "ordinal": 6,
          "snake": "ext_pot_path",
          "type": "Text"
        },
        {
          "capnp": "potentialsPath",
          "default": "",
          "ordinal": 7,
          "snake": "potentials_path",
          "type": "Text"
        }
      ],
      "struct": "PotentialOptions"
    },
    "Process Search": {
      "fields": [
        {
          "capnp": "minimizeFirst",
          "default": True,
          "ordinal": 0,
          "snake": "minimize_first",
          "type": "Bool"
        },
        {
          "capnp": "minimizationOffset",
          "default": 0.2,
          "ordinal": 1,
          "snake": "minimization_offset",
          "type": "Float64"
        }
      ],
      "struct": "ProcessSearchOptions"
    },
    "Structure Comparison": {
      "fields": [
        {
          "capnp": "distanceDifference",
          "default": 0.1,
          "ordinal": 0,
          "snake": "distance_difference",
          "type": "Float64"
        },
        {
          "capnp": "neighborCutoff",
          "default": 3.3,
          "ordinal": 1,
          "snake": "neighbor_cutoff",
          "type": "Float64"
        },
        {
          "capnp": "checkRotation",
          "default": False,
          "ordinal": 2,
          "snake": "check_rotation",
          "type": "Bool"
        },
        {
          "capnp": "indistinguishableAtoms",
          "default": True,
          "ordinal": 3,
          "snake": "indistinguishable_atoms",
          "type": "Bool"
        },
        {
          "capnp": "energyDifference",
          "default": 0.01,
          "ordinal": 4,
          "snake": "energy_difference",
          "type": "Float64"
        },
        {
          "capnp": "removeTranslation",
          "default": True,
          "ordinal": 5,
          "snake": "remove_translation",
          "type": "Bool"
        },
        {
          "capnp": "useCovalent",
          "default": False,
          "ordinal": 6,
          "snake": "use_covalent",
          "type": "Bool"
        },
        {
          "capnp": "covalentScale",
          "default": 1.3,
          "ordinal": 7,
          "snake": "covalent_scale",
          "type": "Float64"
        },
        {
          "capnp": "bruteNeighbors",
          "default": False,
          "ordinal": 8,
          "snake": "brute_neighbors",
          "type": "Bool"
        }
      ],
      "struct": "StructureComparisonOptions"
    }
  },
  "server_yaml_default_overrides": {
    "Main.job": "akmc",
    "Potential.ext_pot_path": "ext_pot"
  },
  "source": "schema/eon_params.capnp"
}
