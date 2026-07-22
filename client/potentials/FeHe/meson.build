# Fortran-only shared module: loaded at runtime via dlopen/LoadLibrary
eon_fehe = shared_module(
    'eon_fehe',
           'feforce.f',
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
                'eon_fehe.def',
            ),
        ] : []
    ),
    install: true,
)

# C++ wrapper (compiled into client, loads Fortran lib at runtime)
fehe_cpp = static_library(
    'fehe_cpp',
    'FeHe.cpp',
    dependencies: _deps,
    cpp_args: _args,
    link_with: _linkto,
    include_directories: _incdirs,
)
