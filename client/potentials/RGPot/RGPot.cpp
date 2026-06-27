/*
** This file is part of eOn.
*/
#include "RGPot.h"
#include "RGPotEngine.h"

RGPot::RGPot(const Parameters &p) : Potential(PotType::RGPot, p) {
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
  impl_ = std::make_unique<RGPotEngine>(opt);
  backend_ = impl_->backend();
}

RGPot::~RGPot() = default;

bool RGPot::engineAvailable() const { return impl_ && impl_->available(); }

void RGPot::force(long N, const double *R, const int *atomicNrs, double *F,
                  double *U, double *variance, const double *box) {
  (void)variance;
  impl_->force(N, R, atomicNrs, F, U, box);
}
