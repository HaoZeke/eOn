Restructure client headers into a real package tree ``include/eon/``
(numpy/fmt/spdlog style). Implementation sources stay under ``client/``.
In-tree and installed form match: ``-I$prefix/include`` and
``#include "eon/Potential.h"`` (and ``eon/fpe_handler.h``,
``eon/potentials/...``). Relative ``../`` includes are gone; headers
install via ``install_subdir`` of that tree (no configure-time header walk).
