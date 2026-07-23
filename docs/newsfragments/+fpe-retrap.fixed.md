FPE continue handler demotes the trapping mask after the first report of each
exception class. Clearing sticky flags alone re-executed the faulting
instruction and re-raised SIGFPE forever (multi-GB identical
`FPE (continuing): division by zero` lines, client stuck at ~100% CPU). Guard
remaining bare divisions and normalizes on the CG / dimer / min-mode saddle
path with SafeMath.
