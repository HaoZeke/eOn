# Fortran-only shared module: loaded at runtime via dlopen/LoadLibrary
eon_edip = shared_module('eon_edip',
           'edipFortran.f90',
           fortran_args : _fargs,
           link_args : _fortran_pot_link_args,
           install : true)

# C++ wrapper (compiled into client, loads Fortran lib at runtime)
edip_cpp = static_library('edip_cpp',
           'EDIP.cpp',
           dependencies : _deps,
           cpp_args : _args,
           link_with : _linkto,
           include_directories: _incdirs)
