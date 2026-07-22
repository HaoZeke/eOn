# Fortran-only shared library: loaded at runtime via dlopen
eon_lenosky = shared_library('eon_lenosky',
           'lenoskyFortran.f90',
           fortran_args : _fargs,
           link_args : _fortran_pot_link_args,
           install : true)

# C++ wrapper (compiled into client, loads Fortran lib at runtime)
lenosky_cpp = static_library('lenosky_cpp',
           'Lenosky.cpp',
           dependencies : _deps,
           cpp_args : _args,
           link_with : _linkto,
           include_directories: _incdirs)
