# cpp_args: _args is required so -DEONMPI / -DWITH_RGPOT reach RGPotEngine.cpp
# (multi-rank PEF force pool). Without it, workers miss collective force under
# mpirun even when the client links MPI.
rgpot_pot = library(
    'rgpot_pot',
    sources: ['RgpotPot.cpp', 'RGPotEngine.cpp'],
    include_directories: _incdirs,
    dependencies: _deps,
    link_with: _linkto,
    cpp_args: _args,
    install: true,
)
