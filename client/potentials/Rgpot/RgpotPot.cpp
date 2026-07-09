/*
** This file is part of eOn.
*/
#include "RgpotPot.h"
#include "RGPotEngine.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {
// Process-lifetime holder for multi-rank rank-0 driver pot (see ClientEON).
std::shared_ptr<RgpotPot> g_rgpot_mpi_driver;
} // namespace

RgpotPot::RgpotPot(const Parameters &p)
    : Potential(PotType::RGPOT, p) {
  RGPotEngineOptions opt;
  const auto &o = p.rgpot_options;
  opt.backend = o.backend;
  opt.basis = o.basis;
  opt.theory = o.theory;
  opt.scf_type = o.scf_type;
  opt.functional = o.functional;
  opt.cutoff_ry = o.cutoff_ry;
  opt.charge = o.charge;
  opt.multiplicity = o.multiplicity;
  opt.engine_path = o.engine_path;
  opt.engine_library = o.engine_library;
  opt.engine_root = o.engine_root;
  opt.title = o.title;
  opt.memory_mb = o.memory_mb;
  opt.scratch_dir = o.scratch_dir;
  opt.input_block = o.input_block;
#ifdef EONMPI
  // Prefer explicit PEF group (mirror LAMMPS MPIClientComm); else WORLD.
  if (p.potential_options.MPIClientComm != MPI_COMM_NULL)
    opt.pe_comm = p.potential_options.MPIClientComm;
  else
    opt.pe_comm = MPI_COMM_WORLD;
#endif

  // Env overrides (CI / benchmarks)
  if (const char *e = std::getenv("RGPOT_BACKEND"))
    opt.backend = e;
  if (const char *e = std::getenv("RGPOT_NWCHEM_BASIS"))
    opt.basis = e;
  if (const char *e = std::getenv("RGPOT_NWCHEM_THEORY"))
    opt.theory = e;
  if (const char *e = std::getenv("RGPOT_NWCHEM_SCF_TYPE"))
    opt.scf_type = e;
  // Engine-path env overrides are backend-scoped: NWCHEMC_LIBRARY must not
  // leak into a cpmdc configure (CPMDPot resolves CPMDC_LIBRARY itself).
  std::string backend_lc = opt.backend;
  std::transform(backend_lc.begin(), backend_lc.end(), backend_lc.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  if (backend_lc.rfind("nwchem", 0) == 0) {
    if (const char *e = std::getenv("NWCHEMC_LIBRARY"))
      opt.engine_path = e;
    else if (const char *e = std::getenv("RGPOT_NWCHEMC_ENGINE"))
      opt.engine_path = e;
  } else if (backend_lc.rfind("cpmd", 0) == 0) {
    if (const char *e = std::getenv("CPMDC_LIBRARY"))
      opt.engine_path = e;
    else if (const char *e = std::getenv("RGPOT_CPMDC_ENGINE"))
      opt.engine_path = e;
  }

  impl_ = std::make_unique<RGPotEngine>(opt);
  backend_ = impl_->backend();
  std::cout << "RgpotPot: in-process rgpot backend=" << backend_
            << " (dlopen libnwchemc/libcpmdc, no potserv RPC)" << std::endl;
#ifdef EONMPI
  const int psz = RGPotEngine::mpiWorldSize();
  if (psz > 1) {
    std::cout << "RgpotPot: multi-rank PEF (EONMPI) world_size=" << psz
              << " rank=" << RGPotEngine::mpiWorldRank() << std::endl;
  }
#endif
}

RgpotPot::~RgpotPot() = default;

bool RgpotPot::engineAvailable() const { return impl_ && impl_->available(); }

void RgpotPot::force(long N, const double *R, const int *atomicNrs, double *F,
                     double *U, double *variance, const double *box) {
  (void)variance;
  impl_->force(N, R, atomicNrs, F, U, box);
}

void RgpotPot::mpiWorkerLoop() const {
  if (!impl_)
    throw std::runtime_error("RgpotPot::mpiWorkerLoop: no engine");
  impl_->mpiWorkerLoop();
}

void RgpotPot::mpiSignalStop() const {
  if (impl_)
    impl_->mpiSignalStop();
}

int RgpotPot::mpiWorldSize() { return RGPotEngine::mpiWorldSize(); }
int RgpotPot::mpiWorldRank() { return RGPotEngine::mpiWorldRank(); }

void RgpotPot::setMpiDriverInstance(std::shared_ptr<RgpotPot> pot) {
  g_rgpot_mpi_driver = std::move(pot);
}

std::shared_ptr<RgpotPot> RgpotPot::mpiDriverInstance() {
  return g_rgpot_mpi_driver;
}

void RgpotPot::clearMpiDriverInstance() { g_rgpot_mpi_driver.reset(); }
