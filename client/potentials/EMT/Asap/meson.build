# Headers live under include/eon/potentials/EMT/Asap/; short includes
# (Atoms.h, EMT.h, ...) need that directory on the include path.
_asap_hdr_inc = include_directories(
    join_paths(
        meson.project_source_root(),
        'include',
        'eon',
        'potentials',
        'EMT',
        'Asap',
    ),
)
asap = library(
    'asap',
    [
        'NeighborList.cpp',
        'GhostAtoms.cpp',
        'CellList.cpp',
        'Vec.cpp',
        'Atoms.cpp',
        'Exception.cpp',
        'EMTRasmussenParameterProvider.cpp',
        'EMTDefaultParameterProvider.cpp',
        'EMT.cpp',
        'stacktrace.c',
        'SuperCell.cpp',
    ],
    dependencies: _deps,
    link_with: _linkto,
    cpp_args: _args,
    include_directories: _incdirs + [_asap_hdr_inc],
    install: true,
)
