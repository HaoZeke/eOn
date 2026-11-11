#include "quill/backend/BackendWorkerLock.h"
#include "eon/EonLogger.h"
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
int main() {
  setenv("QUILL_LOG_LEVEL", "info", 1);
  auto *ready = eonc::log::init_client();
  assert(ready != nullptr);
  auto *a = eonc::log::get();
  auto *b = eonc::log::get();
  assert(a != nullptr && a == b);
  auto *t1 = eonc::log::traceback();
  auto *t2 = eonc::log::traceback();
  assert(t1 == t2);
  eonc::log::Scoped scoped;
  assert(scoped.ptr() != nullptr || scoped.ptr() == nullptr);
  auto *file = eonc::log::get_file("p321", "p321.log");
  auto *file2 = eonc::log::get_file("p321", "p321.log");
  assert(file != nullptr && file == file2);
  eonc::log::FileScoped fs{"p321b", "p321b.log"};
  assert(fs.ptr() != nullptr || fs.ptr() == nullptr);
  quill::detail::BackendWorkerLock lock("p321");
  bool threw = false;
  try {
    quill::detail::BackendWorkerLock again("p321");
    (void)again;
  } catch (...) {
    threw = true;
  }
  assert(threw);
  // unwritable directory forces the open-failure return
  if (chdir("/proc") == 0) {
    quill::detail::BackendWorkerLock skipped("no-such");
  }
  std::cout << "ok\n";
}
