Client C++ sources use logical includes from the `client/` root
(`#include "Potential.h"`, `#include "potentials/LAMMPS/LAMMPSPot.h"`)
instead of relative `../` paths. Headers install under the FHS path
`$prefix/include/eonclient/` with the same layout.
