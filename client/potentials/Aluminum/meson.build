# Fortran-only shared module: loaded at runtime via dlopen/LoadLibrary
eon_aluminum = shared_module(
    'eon_aluminum',
           'fofrhoDblexp.f',
           'dfrhoDblexp.f',
           'potinit.f',
           'forces.f',
           'gagafeDblexp.f',
           'embedenergy.f',
           'sumembforce.f',
    fortran_args: _fargs,
    link_args: _fortran_pot_link_args
    + (
        is_windows ? [
            '/DEF:'
            + join_paths(
                meson.project_source_root(),
                'client',
                'potentials',
                'win_exports',
                'eon_aluminum.def',
            ),
        ] : []
    ),
    install: true,
)

# C++ wrapper (compiled into client, loads Fortran lib at runtime)
aluminum_cpp = static_library(
    'aluminum_cpp',
    'Aluminum.cpp',
    dependencies: _deps,
    cpp_args: _args,
    link_with: _linkto,
    include_directories: _incdirs,
)
