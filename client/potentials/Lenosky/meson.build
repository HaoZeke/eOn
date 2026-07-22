# Fortran-only shared module: loaded at runtime via dlopen/LoadLibrary
eon_lenosky = shared_module(
    'eon_lenosky',
           'lenoskyFortran.f90',
    fortran_args: _fargs,
    link_args: _fortran_pot_link_args
    + (
        is_windows ? [
            '-Wl,/DEF:'
            + join_paths(
                meson.project_source_root(),
                'client',
                'potentials',
                'win_exports',
                'eon_lenosky.def',
            ),
        ] : []
    ),
    install: true,
)

# C++ wrapper (compiled into client, loads Fortran lib at runtime)
lenosky_cpp = static_library(
    'lenosky_cpp',
    'Lenosky.cpp',
    dependencies: _deps,
    cpp_args: _args,
    link_with: _linkto,
    include_directories: _incdirs,
)
