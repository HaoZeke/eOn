LAMMPS worker process demotes floating-point traps after fork, and
`forceLocal` holds an `FPEHandler::eat_fpe` scope. Client `enableFPE()` was
inherited by the forked PairEAM worker; a single FE_DIVBYZERO inside
`LAMMPS_NS::PairEAM::compute` then re-stormed under the broken continue
handler (multi-GB identical stderr lines, child at ~100% CPU). Matches the
ASE/Metatomic external-pot contract.
