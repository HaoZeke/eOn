# Built AFTER eclib: torch+metatomic plugin for dlopen. Host never links this.
mta_libdir = meson.current_source_dir()
metatomic_pot = library(
    'metatomic_pot',
    [
        'MetatomicPotential.cpp',
        'MetatomicCAbi.cpp',
    ],
    dependencies: mta_plugin_deps,
    cpp_args: _args + ['-DEON_MTA_BUILD'],
    # Link full host libs so the plugin is self-contained (destroy may still
    # need RTLD_GLOBAL binding; see MetatomicLoader).
    link_with: [eclib, eoncbase],
    include_directories: _incdirs,
    build_rpath: ':'.join(_runtime_rpath),
    install_rpath: ':'.join(_runtime_rpath),
    install: true,
    gnu_symbol_visibility: 'default',
)
