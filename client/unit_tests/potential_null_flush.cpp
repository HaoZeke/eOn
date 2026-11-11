/*
** Copyright (c) eOn developers
*/

// A bad QUILL_LOG_LEVEL makes get() return null.
// The unknown-potential arm used to call flush_log on that null.
// Removing the null check aborts this process, so the test fails.

#include "eon/Parameters.h"
#include "eon/Potential.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

int main() {
  ::setenv("QUILL_LOG_LEVEL", "not-a-level", 1);
  eonc::Parameters params;
  try {
    auto pot = eonc::makePotential(eonc::PotType::UNKNOWN, params);
    if (pot) {
      std::cerr << "unknown potential was constructed\n";
      return 1;
    }
  } catch (const std::runtime_error &err) {
    if (std::string(err.what()) == "Terminating") {
      std::cout << "null_flush caught\n";
      return 0;
    }
    std::cerr << err.what() << "\n";
    return 2;
  } catch (...) {
    std::cerr << "unexpected exception\n";
    return 3;
  }
  return 1;
}
