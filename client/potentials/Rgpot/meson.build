# cpp_args: _args is required so -DWITH_RGPOT_MPI / -DWITH_RGPOT reach
# RGPotEngine.cpp (multi-rank force pool). Without it, workers throw
# "built without WITH_RGPOT_MPI" under mpirun even when the client links MPI.
rgpot_pot = library(
    'rgpot_pot',
    sources: ['RgpotPot.cpp', 'RGPotEngine.cpp'],
    include_directories: _incdirs,
    dependencies: _deps,
    link_with: _linkto,
    cpp_args: _args,
    install: true,
)
