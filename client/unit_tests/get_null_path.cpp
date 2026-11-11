/*
** Copyright (c) eOn developers
*/

// A bad QUILL_LOG_LEVEL makes create_or_get_logger throw.
// get() catches that throw and returns a null pointer.
// Removing the catch aborts this process, so the test fails.

#include "eon/EonLogger.h"

#include <cstdlib>
#include <iostream>
#include <string>

int main(int argc, char **argv) {
  const std::string mode = argc > 1 ? argv[1] : "bad-level";
  if (mode == "bad-level") {
    ::setenv("QUILL_LOG_LEVEL", "not-a-level", 1);
    quill::Logger *got = nullptr;
    try {
      got = eonc::log::get();
    } catch (...) {
      std::cerr << "exception escaped get()\n";
      return 2;
    }
    if (got != nullptr) {
      std::cerr << "get() returned a logger\n";
      return 1;
    }
    std::cout << "bad_level get 0\n";
    return 0;
  }
  std::cerr << "unknown mode\n";
  return 3;
}
