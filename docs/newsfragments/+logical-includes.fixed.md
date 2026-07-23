Client C++ uses package-prefixed includes (`#include "eon/Potential.h"`,
`#include "eon/fpe_handler.h"`, `#include "eon/potentials/LAMMPS/LAMMPSPot.h"`)
with in-tree `-I$src/include` (`include/eon` -> `client`) and FHS install under
`$prefix/include/eon/` (same form as `xtb/xtb.h` with `-I$prefix/include`).
Relative `../` includes are gone.
