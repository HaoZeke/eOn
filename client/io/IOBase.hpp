/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#pragma once
#include "client/Matter.h"

#include <fstream>
#include <string>

namespace eonc::io {

// Base class for common IO operations
class IOBase {
public:
  bool fileExists(const std::string &);
  void ensureFileOpen(std::ofstream &, const std::string &, bool);
  void ensureFileOpen(std::ifstream &, const std::string &);
};

// Base class for writing
template <typename Derived> class WriteBase : public IOBase {
public:
  bool write(const Matter &matter, const std::string &filename,
             bool append = false) {
    std::ofstream file;
    ensureFileOpen(file, filename, append);

    if (!file.is_open()) {
      return false;
    }

    bool success = static_cast<Derived *>(this)->writeImpl(matter, file);
    file.close();
    return success;
  }
};

// Base class for reading
template <typename Derived> class ReadBase : public IOBase {
public:
  bool read(Matter &matter, const std::string &filename) {
    std::ifstream file;
    ensureFileOpen(file, filename);

    if (!file.is_open()) {
      return false;
    }

    bool success = static_cast<Derived *>(this)->readImpl(matter, file);
    file.close();
    return success;
  }
};

} // namespace eonc::io
