/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
*/

// Four public entries to get(). Each prints 1 when the pointer is set.
// A zero would mean the log-null arm is reachable on that entry.

#include "eon/EonLogger.h"

#include <iostream>
#include <string>

int main(int argc, char **argv) {
  std::string mode = argc > 1 ? argv[1] : "fresh";
  if (mode == "init") {
    auto *created = eonc::log::init_client();
    auto *got = eonc::log::get();
    std::cout << "init " << (created ? 1 : 0) << " get " << (got ? 1 : 0)
              << "\n";
    return got ? 0 : 2;
  }
  if (mode == "stop") {
    auto *created = eonc::log::init_client();
    quill::Backend::stop();
    auto *got = eonc::log::get();
    std::cout << "after_stop init " << (created ? 1 : 0) << " get "
              << (got ? 1 : 0) << "\n";
    return got ? 0 : 3;
  }
  if (mode == "stop-first") {
    quill::Backend::stop();
    auto *got = eonc::log::get();
    std::cout << "stop_first get " << (got ? 1 : 0) << "\n";
    return got ? 0 : 4;
  }
  auto *got = eonc::log::get();
  std::cout << "fresh get " << (got ? 1 : 0) << "\n";
  EONC_LOG_DEBUG("fresh line");
  return got ? 0 : 1;
}
