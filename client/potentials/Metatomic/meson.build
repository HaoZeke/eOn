# Fat pot library: linked into eonclient via potentials (potential = Metatomic).
metatomic_pot = library(
    'metatomic_pot',
    'MetatomicPotential.cpp',
    dependencies: _deps,
    cpp_args: _args,
    link_with: _linkto,
    include_directories: _incdirs,
    install: true,
)

# Optional engine for RGPOT-metatomic dlopen (same MetatomicPotential force path).
# Built only when eclib already exists (after eclib in client/meson.build).
if is_variable('eclib')
    metatomic_engine = library(
        'metatomic_engine',
        [
            'MetatomicPotential.cpp',
            'MetatomicEngineAbi.cpp',
        ],
        dependencies: mta_plugin_deps,
        cpp_args: _args + ['-DRGPOT_MTA_ENGINE_BUILD'],
        link_with: [eclib],
        include_directories: _incdirs,
        build_rpath: ':'.join(_runtime_rpath),
        install_rpath: ':'.join(_runtime_rpath),
        install: true,
        gnu_symbol_visibility: 'default',
    )
endif
