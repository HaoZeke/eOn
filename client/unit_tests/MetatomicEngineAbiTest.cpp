#include "catch2/catch_amalgamated.hpp"
#include <dlfcn.h>
#include <string>

// Structural + optional live check: libmetatomic_engine exports rgpot_mta_* ABI.
TEST_CASE("metatomic engine C ABI symbols present when engine is loadable",
          "[metatomic][rgpot][engine]") {
  const char *env = std::getenv("RGPOT_METATOMIC_ENGINE");
  std::string path = env && *env ? env : "libmetatomic_engine.so";
  void *h = dlopen(path.c_str(), RTLD_NOW | RTLD_GLOBAL);
  if (!h) {
    // Also try build-tree relative path from test env
    const char *pp = std::getenv("EON_POTENTIALS_PATH");
    if (pp) {
      std::string d(pp);
      auto pos = d.find(':');
      if (pos != std::string::npos)
        d = d.substr(0, pos);
      if (!d.empty() && d.back() != '/')
        d += '/';
      h = dlopen((d + "libmetatomic_engine.so").c_str(), RTLD_NOW | RTLD_GLOBAL);
    }
  }
  if (!h) {
    WARN("libmetatomic_engine.so not on path; skip live ABI check");
    return;
  }
  REQUIRE(dlsym(h, "rgpot_mta_abi_version") != nullptr);
  REQUIRE(dlsym(h, "rgpot_mta_create") != nullptr);
  REQUIRE(dlsym(h, "rgpot_mta_destroy") != nullptr);
  REQUIRE(dlsym(h, "rgpot_mta_force") != nullptr);
  auto abi = reinterpret_cast<int (*)()>(dlsym(h, "rgpot_mta_abi_version"));
  REQUIRE(abi != nullptr);
  REQUIRE(abi() == 1);
  dlclose(h);
}
