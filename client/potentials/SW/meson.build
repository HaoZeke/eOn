# Fortran-only shared module: loaded at runtime via dlopen/LoadLibrary
eon_sw = shared_module(
    'eon_sw',
           'SWFortran.f90',
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
                'eon_sw.def',
            ),
        ] : []
    ),
    install: true,
)

# C++ wrapper (compiled into client, loads Fortran lib at runtime)
sw_cpp = static_library(
    'sw_cpp',
    'SW.cpp',
    dependencies: _deps,
    cpp_args: _args,
    link_with: _linkto,
    include_directories: _incdirs,
)
