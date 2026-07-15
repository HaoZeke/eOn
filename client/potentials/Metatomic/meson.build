# Built AFTER eclib (see client/meson.build): full link of torch + eonclib so the
# shared plugin is self-sufficient for dlopen. Host never links this target.
mta_libdir = meson.current_source_dir()
metatomic_pot = library(
    'metatomic_pot',
    [
        'MetatomicPotential.cpp',
        'MetatomicCAbi.cpp',
    ],
    dependencies: mta_plugin_deps,
    cpp_args: _args + ['-DEON_MTA_BUILD'],
    # Host already provides Potential/PotRegistry via libeonclib; leave those
    # unresolved so dlopen(RTLD_GLOBAL) binds to the process-wide symbols.
    link_with: [eoncbase],
    link_args: ['-Wl,--allow-shlib-undefined'],
    include_directories: _incdirs,
    build_rpath: ':'.join(_runtime_rpath),
    install_rpath: ':'.join(_runtime_rpath),
    install: true,
    gnu_symbol_visibility: 'default',
)
