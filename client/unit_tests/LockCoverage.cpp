// Copyright (c) amsel developers
#include <cassert>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <new>
#include <string>
#include <stdexcept>
#include <sys/file.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <vector>

static std::vector<int> open_codes;
static std::vector<int> flock_codes;
static size_t open_i, flock_i;
static int live_fd = -1;
static int throw_after = -1;
static int news = 0;

void *operator new(std::size_t n) {
  if (throw_after >= 0 && news++ == throw_after) {
    throw std::bad_alloc();
  }
  return std::malloc(n);
}
void operator delete(void *p) noexcept { std::free(p); }
void operator delete(void *p, std::size_t) noexcept { std::free(p); }

static int sys_open(const char *path, int flags, int mode) {
  return syscall(SYS_openat, AT_FDCWD, path, flags, mode);
}

extern "C" int p322_open(const char *path, int flags, ...) {
  int code = open_codes.at(open_i++);
  if (code == -3) {
    throw std::runtime_error("open");
  }
  if (code == -4) {
    errno = EINTR;
    return -1;
  }
  if (code < 0) {
    errno = EIO;
    return -1;
  }
  live_fd = sys_open(path, flags, 0644);
  assert(live_fd >= 0);
  return live_fd;
}
extern "C" int p322_flock(int fd, int op) {
  int code = flock_codes.at(flock_i++);
  if (code == -3) {
    throw std::runtime_error("flock");
  }
  if (code == -4) {
    errno = EINTR;
    return -1;
  }
  if (code == -2) {
    errno = EWOULDBLOCK;
    return -1;
  }
  if (code < 0) {
    errno = EIO;
    return -1;
  }
  return 0;
}
extern "C" int p322_close(int fd) {
  if (fd == live_fd && live_fd >= 0) {
    int rc = syscall(SYS_close, fd);
    live_fd = -1;
    return rc;
  }
  return 0;
}
extern "C" int p322_unlink(const char *path) {
  return syscall(SYS_unlinkat, AT_FDCWD, path, 0);
}
inline void p322_sleep(unsigned long) {}

#define open p322_open
#define flock p322_flock
#define close p322_close
#define unlink p322_unlink
#define sleep_for_ns p322_sleep
#include "quill/backend/BackendWorkerLock.h"
#undef open
#undef flock
#undef close
#undef unlink
#undef sleep_for_ns

static void run(const std::vector<int> &opens, const std::vector<int> &flocks,
                bool expect_throw, const char *label) {
  open_codes = opens;
  flock_codes = flocks;
  open_i = flock_i = 0;
  bool threw = false;
  try {
    quill::detail::BackendWorkerLock lock(label);
  } catch (const std::exception &) {
    threw = true;
  }
  assert(threw == expect_throw);
  assert(open_i == open_codes.size());
}

int main() {
  // success, destructor releases
  run({0}, {0}, false, "ok");
  // EINTR then success on open; flock success
  run({-4, 0}, {0}, false, "eintr-open");
  // open fails with EIO three times and the constructor returns
  run({-1, -1, -1}, {}, false, "open-fail");
  // flock EINTR then success
  run({0}, {-4, 0}, false, "eintr-flock");
  // duplicate lock throws
  run({0}, {-2}, true, "dup");
  // flock EIO retries then gives up without throwing
  run({0}, {-1, -1, -1}, false, "flock-io");
  // the open and flock calls themselves throw
  run({-3}, {}, true, "open-throws");
  run({0}, {-3}, true, "flock-throws");
  open_codes = {0};
  flock_codes = {0};
  for (int n = 0; n < 6; ++n) {
    news = 0;
    throw_after = n;
    open_i = flock_i = 0;
    try {
      quill::detail::BackendWorkerLock lock("alloc");
    } catch (const std::bad_alloc &) {
    } catch (const std::exception &) {
    }
  }
  // Throw while the duplicate-lock error string is built.
  open_codes = {0};
  flock_codes = {-2};
  for (int n = 0; n < 12; ++n) {
    news = 0;
    throw_after = n;
    open_i = flock_i = 0;
    try {
      quill::detail::BackendWorkerLock lock("dupstr");
    } catch (const std::bad_alloc &) {
    } catch (const std::exception &) {
    }
  }
  throw_after = -1;
  std::cout << "lock-ok\n";
}
