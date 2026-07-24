``[Hessian] phva_atoms`` replaces ``atom_list`` for the PHVA mobile set
(default still ``All``). The C++ client still accepts the legacy
``atom_list`` key when ``phva_atoms`` is absent. Lanczos and Davidson
gain the same ``phva_atoms`` key in INI, schema, and ``config.yaml``.
