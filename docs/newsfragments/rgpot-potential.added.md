Added an ``RGPot`` client potential that drives the rgpot ``NWChemPot`` and
``CPMDPot`` backends in-process (via the rgpot C ABI / dlopen engines), as a
faster alternative to the socket-based ``SocketNWChem`` potential. Enable with
Meson ``-Dwith_rgpot=true`` (or CMake ``-DWITH_RGPOT``); select the backend with
the ``[RGPot]`` ``backend`` option (``nwchemc`` or ``cpmdc``). Includes a
CPMD-BLYP point-energy example under ``examples/rgpot-cpmd-blyp/``.
