#pragma once

#include <memory>
#include <string>

#ifdef EONMPI
#include <mpi.h>
#endif

struct RGPotEngineOptions {
  std::string backend;
  std::string basis{"sto-3g"};
  std::string theory{"scf"};
  std::string scf_type{"rhf"};
  std::string functional{"BLYP"};
  double cutoff_ry{70.0};
  int charge{0};
  int multiplicity{1};
  std::string engine_path;
  std::string engine_library;
  std::string engine_root;
  std::string title;
  int memory_mb{0};
  std::string scratch_dir;
  std::string input_block; // optional NWChem inputBlocks text
#ifdef EONMPI
  /// PEF communicator (host-owned MPI). MPI_COMM_NULL → MPI_COMM_WORLD.
  MPI_Comm pe_comm{MPI_COMM_NULL};
#endif
};

/** Opaque rgpot-backed engine (nwchemc / cpmdc). No eOn Potential.h here. */
class RGPotEngine {
public:
  explicit RGPotEngine(const RGPotEngineOptions &opt);
  ~RGPotEngine();
  RGPotEngine(const RGPotEngine &) = delete;
  RGPotEngine &operator=(const RGPotEngine &) = delete;

  [[nodiscard]] const std::string &backend() const noexcept { return backend_; }
  [[nodiscard]] bool available() const;

  /**
   * Rank 0 of the PEF communicator owns N/R/Z/box from eOn. When the PEF
   * communicator has size > 1, force() is collective: broadcast geometry,
   * every rank evaluates the engine, rank 0 keeps U/F for eOn.
   */
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, const double *box) const;

  /**
   * Multi-rank PEF pool: non-zero ranks block here and join the same
   * collectives as force() until rank 0 broadcasts N<=0 (shutdown).
   */
  void mpiWorkerLoop() const;

  /** World size if MPI is active, else 1. */
  static int mpiWorldSize();
  static int mpiWorldRank();

  /// Rank 0 only: broadcast N=0 so workers leave mpiWorkerLoop().
  void mpiSignalStop() const;

#ifdef EONMPI
  [[nodiscard]] MPI_Comm peComm() const noexcept { return pe_comm_; }
#endif

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
  std::string backend_;
#ifdef EONMPI
  MPI_Comm pe_comm_{MPI_COMM_WORLD};
#endif

  void forceLocal(long N, const double *R, const int *atomicNrs, double *F,
                  double *U, const double *box) const;
};
