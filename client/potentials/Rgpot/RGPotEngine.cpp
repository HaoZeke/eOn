// Isolated TU: rgpot only (no eOn Potential.h) — avoids Cap'n Proto Potential
// clash.
#include "RGPotEngine.h"

#include <array>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <vector>

#include <capnp/message.h>
#include <capnp/serialize.h>

#include "rgpot/CPMDPot/CPMDPot.hpp"
#include "rgpot/NWChemPot/NWChemPot.hpp"
#include "rgpot/rpc/Potentials.capnp.h"

namespace {

std::string to_lower(std::string s) {
  for (char &c : s)
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

std::array<std::array<double, 3>, 3> box_from_row_major(const double *box) {
  std::array<std::array<double, 3>, 3> out{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      out[static_cast<size_t>(i)][static_cast<size_t>(j)] = box[i * 3 + j];
  return out;
}

bool looks_like_dft_xc(const std::string &s) {
  if (s.empty())
    return false;
  static const char *k[] = {"b3lyp", "blyp", "pbe",    "pw91",    "bp86",
                            "hcth",  "ft97", "hfexch", "xperpbe", nullptr};
  for (int i = 0; k[i]; ++i) {
    if (s.size() >= std::char_traits<char>::length(k[i]) &&
        s.compare(0, std::char_traits<char>::length(k[i]), k[i]) == 0)
      return true;
  }
  return false;
}

// Serialized Cap'n Proto message (standard segment framing, word-aligned) as
// written by `capnp convert text:binary` or messageToFlatArray.
std::vector<::capnp::word> read_params_file(const std::string &path) {
  std::ifstream in(path, std::ios::binary | std::ios::ate);
  if (!in)
    throw std::runtime_error("RGPOT: cannot open params_path: " + path);
  const std::streamsize bytes = in.tellg();
  if (bytes <= 0 ||
      (static_cast<size_t>(bytes) % sizeof(::capnp::word)) != 0)
    throw std::runtime_error("RGPOT: params_path is not a capnp flat message: " +
                             path);
  std::vector<::capnp::word> words(static_cast<size_t>(bytes) /
                                   sizeof(::capnp::word));
  in.seekg(0);
  in.read(reinterpret_cast<char *>(words.data()), bytes);
  if (!in)
    throw std::runtime_error("RGPOT: short read on params_path: " + path);
  return words;
}

} // namespace

struct RGPotEngine::Impl {
  enum class Backend { Nwchemc, Cpmdc };
  Backend backend{Backend::Nwchemc};
  std::unique_ptr<rgpot::NWChemPot> nwchem;
  std::unique_ptr<rgpot::CPMDPot> cpmd;
};

RGPotEngine::RGPotEngine(const RGPotEngineOptions &opt)
    : impl_(std::make_unique<Impl>()) {
  backend_ = to_lower(opt.backend);
  if (backend_ == "nwchem" || backend_ == "nwchemc" ||
      backend_ == "nwchempot") {
    backend_ = "nwchemc";
    impl_->backend = Impl::Backend::Nwchemc;
    ::capnp::MallocMessageBuilder msg;
    auto params = msg.initRoot<::NWChemParams>();
    params.setBasis(opt.basis);
    params.setTheory(opt.theory);
    params.setScfType(opt.scf_type);
    params.setCharge(opt.charge);
    params.setMultiplicity(opt.multiplicity);
    if (!opt.engine_path.empty())
      params.setEnginePath(opt.engine_path);
    else if (!opt.engine_library.empty())
      params.setEnginePath(opt.engine_library);
    if (!opt.engine_root.empty())
      params.setNwchemRoot(opt.engine_root);
    if (!opt.title.empty())
      params.setTitle(opt.title);
    if (opt.memory_mb > 0)
      params.setMemoryMb(static_cast<uint32_t>(opt.memory_mb));
    if (!opt.scratch_dir.empty())
      params.setScratchDir(opt.scratch_dir);
    // DFT XC as inputBlocks when theory=dft and scfType is an XC label
    std::string block = opt.input_block;
    if (block.empty()) {
      if (const char *env = std::getenv("RGPOT_NWCHEM_INPUT_BLOCK"))
        block = env;
    }
    if (block.empty() && (opt.theory == "dft" || opt.theory == "DFT") &&
        looks_like_dft_xc(opt.scf_type)) {
      block = "dft\n  xc " + opt.scf_type + "\n  mult " +
              std::to_string(opt.multiplicity) + "\nend";
    } else if (block.empty() && looks_like_dft_xc(opt.theory)) {
      block = "dft\n  xc " + opt.theory + "\n  mult " +
              std::to_string(opt.multiplicity) + "\nend";
    }
    if (!block.empty()) {
      auto blocks = params.initInputBlocks(1);
      blocks.set(0, block);
    }
    impl_->nwchem = std::make_unique<rgpot::NWChemPot>(params.asReader());
    if (!impl_->nwchem->available())
      throw std::runtime_error(
          "RGPOT(nwchemc): engine not available (set NWCHEMC_LIBRARY / "
          "RGPOT_NWCHEMC_ENGINE or [RgpotPot] engine_path)");
  } else if (backend_ == "cpmd" || backend_ == "cpmdc" ||
             backend_ == "cpmdpot") {
    backend_ = "cpmdc";
    impl_->backend = Impl::Backend::Cpmdc;
    ::capnp::MallocMessageBuilder msg;
    ::CPMDParams::Builder params = msg.initRoot<::CPMDParams>();
    if (!opt.params_path.empty()) {
      // Full CPMDParams message from disk (functional, cutoff, inputSections
      // with method directives and pseudopotentials). INI scalars below only
      // fill engine location / scratch, never the method sections.
      const auto words = read_params_file(opt.params_path);
      ::capnp::FlatArrayMessageReader reader(
          kj::arrayPtr(words.data(), words.size()));
      msg.setRoot(reader.getRoot<::CPMDParams>());
      params = msg.getRoot<::CPMDParams>();
    } else {
      params.setFunctional(opt.functional);
      params.setCutOffRy(opt.cutoff_ry);
      params.setCharge(opt.charge);
      params.setMultiplicity(opt.multiplicity);
      if (!opt.title.empty())
        params.setTitle(opt.title);
      if (opt.memory_mb > 0)
        params.setMemoryMb(static_cast<uint32_t>(opt.memory_mb));
    }
    if (!opt.engine_path.empty())
      params.setEnginePath(opt.engine_path);
    else if (!opt.engine_library.empty())
      params.setEnginePath(opt.engine_library);
    if (!opt.engine_root.empty())
      params.setCpmdRoot(opt.engine_root);
    if (!opt.scratch_dir.empty())
      params.setScratchDir(opt.scratch_dir);
    impl_->cpmd = std::make_unique<rgpot::CPMDPot>(params.asReader());
    if (!impl_->cpmd->available())
      throw std::runtime_error(
          "RGPOT(cpmdc): engine not available (set CPMDC_LIBRARY / "
          "RGPOT_CPMDC_ENGINE or [RgpotPot] engine_path)");
  } else {
    throw std::runtime_error("RGPOT: unknown backend '" + opt.backend +
                             "' (expected nwchemc or cpmdc)");
  }
}

RGPotEngine::~RGPotEngine() = default;

bool RGPotEngine::available() const {
  if (!impl_)
    return false;
  if (impl_->backend == Impl::Backend::Nwchemc && impl_->nwchem)
    return impl_->nwchem->available();
  if (impl_->backend == Impl::Backend::Cpmdc && impl_->cpmd)
    return impl_->cpmd->available();
  return false;
}

void RGPotEngine::force(long N, const double *R, const int *atomicNrs,
                        double *F, double *U, const double *box) const {
  if (N <= 0)
    throw std::runtime_error("RGPotEngine::force called with N <= 0");

  AtomMatrix positions(static_cast<int>(N), 3);
  for (long i = 0; i < N; ++i) {
    const int ii = static_cast<int>(i);
    positions(ii, 0) = R[3 * i + 0];
    positions(ii, 1) = R[3 * i + 1];
    positions(ii, 2) = R[3 * i + 2];
  }
  std::vector<int> atmtypes(atomicNrs, atomicNrs + N);
  const auto cell = box_from_row_major(box);

  std::pair<double, AtomMatrix> result;
  if (impl_->backend == Impl::Backend::Nwchemc)
    result = (*impl_->nwchem)(positions, atmtypes, cell);
  else
    result = (*impl_->cpmd)(positions, atmtypes, cell);

  *U = result.first;
  const auto &forces = result.second;
  for (long i = 0; i < N; ++i) {
    const int ii = static_cast<int>(i);
    F[3 * i + 0] = forces(ii, 0);
    F[3 * i + 1] = forces(ii, 1);
    F[3 * i + 2] = forces(ii, 2);
  }
}
