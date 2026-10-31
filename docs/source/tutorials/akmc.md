---
myst:
  html_meta:
    "description": "Tutorial: Adaptive Kinetic Monte Carlo (aKMC) on a Pt heptamer island on Pt(111)."
    "keywords": "eOn aKMC tutorial, Kinetic Monte Carlo, Pt(111), saddle search, simulation"
---

# AKMC Tutorial

This tutorial runs adaptive kinetic Monte Carlo (aKMC) on a seven-atom platinum island, a heptamer, on a Pt(111) surface.
It stops after the first kinetic Monte Carlo (KMC) step.
The potential is the built-in `morse_pt`, so an eOn installation is all it needs (see the [installation instructions](project:../install/index.md)).
The run shown here took a little over three minutes on an eight-core machine.
The method and its references are on the [aKMC user guide page](project:../user_guide/akmc.md).

## How a run proceeds

```{figure} ../fig/akmc/akmc_cycle.png
---
alt: One server call registers results and updates the process table. Below the confidence target it makes new searches. At the target it cancels queued work and takes a KMC step.
class: full-width
align: center
---
One call of the server.
Results are read one call after the searches run.
```

Each call of `python -m eon.server` first registers the results of the searches made by the previous call.
It adds new processes to the table of processes that leave the current state and updates the confidence that the table is complete.
While that confidence is below the target, the call makes new searches.
With the local communicator used here, the same call starts `eonclient` for each search and waits for it to finish.
When the confidence reaches the target, the call cancels queued work and takes a KMC step: it picks one process with a probability proportional to its rate and moves to that process's product state.
It does not make new searches on that call.
Each search displaces atoms near the least coordinated atom, climbs to a saddle point with the dimer method, minimizes on both sides of the saddle and computes the rate prefactor.

## Set up the run

Copy the example out of the source tree and work in the copy:

```{code-block} bash
cp -r examples/akmc-pt akmc-pt-run
cd akmc-pt-run
```

`pos.con` holds the 343 atoms of the slab and the island.
`config.ini` is:

```{literalinclude} ../../../examples/akmc-pt/config.ini
:language: ini
```

`random_seed = 1024` makes the server draw the same random numbers on every run, so a run on the same build repeats the numbers on this page.
Another build or machine can give other numbers.
`confidence = 0.95` is the target for the KMC step.
`number_of_CPUs = 2` and `num_jobs = 2` make each call run two searches side by side.

## The first two calls

```{code-block} console
$ python -m eon.server
Eon version 3.6.0
State list path does not exist; Creating: .//states/
Registering results
Processed 0 results
Queue contains 0 searches
Making 2 process searches
Job finished: jobs/scratch/0_1
Job finished: jobs/scratch/0_0
Created 2 searches
Currently in state 0 with confidence 0.000000
```

The first call has nothing to register.
It makes two searches and runs them in `jobs/scratch/`.
The order of the `Job finished` lines can change from run to run.
The second call registers those two results and makes two more.
Neither of the first two searches found a process.
The third call prints `Found new lowest barrier 0.601061 for state 0`, the first process.

Each search adds one row to `states/0/search_results.txt`.
These are its first rows:

```{code-block} text
    wuid       type    barrier   max-dist    sad-fcs   mins-fcs   pref-fcs    result
-------------------------------------------------------------------------------------
       0     random    0.00000    5.39201        701        157          0    Not Connected
       1     random    0.00000    2.82617        423        259          0    Not Connected
       3     random    0.00000    2.07997        379        251          0    Not Connected
       2     random    0.60106    0.94440        305        187       1404    good-1905950264122117683
       5     random    1.97791    3.68020        326        303       1278    barrier > max_thermal_window
       4     random    0.00000    0.90982        154          1          0    repeat-1905950264122117683
```

`good-<id>` is a new process and `repeat-<id>` is a process found before.
`Not Connected` means that the two minimizations from the saddle did not give the current state on one side and a different state on the other.
`barrier > max_thermal_window` marks a saddle so far above the lowest barrier that it is dropped.
The `-fcs` columns count force calls for the saddle search, the two minimizations and the prefactor.

## Run to the first KMC step

Call the server until it reports a KMC step in `akmc.log`:

```{code-block} bash
until grep -q "KMC step from" akmc.log; do python -m eon.server; done
```

```{figure} ../fig/akmc/akmc_pt_confidence.png
---
alt: Confidence in state 0 against the number of search results registered
align: center
---
Confidence printed after each call.
It falls back to zero each time a new process is found.
```

With the default `confidence_scheme = old`, the confidence is {math}`1 - 1/N`.
{math}`N` is the number of searches that found a process already in the rate table, counted since the last new process.
The rate table holds the processes whose barriers are within a thermal window of the lowest barrier.
A new process sets {math}`N` back to zero, which is why the curve falls.
A target of 0.95 needs 20 repeats with no new process between them.
In this run the confidence first rose above zero on call 13.
The KMC step came on call 74, after 146 search results:

```{code-block} console
$ python -m eon.server
Eon version 3.6.0
Registering results
Processed 2 results
Cancelled 0 workunits from state 0
KMC step from state 0 through process 1905950264122117683 to state 1
Currently in state 1 with confidence 0.000000
```

## The process table

`states/0/processtable` lists the 8 processes found in state 0:

```{code-block} text
 proc #    saddle energy   prefactor   product   product energy product prefactor  barrier         rate repeats
1905950264122117683      -1775.19009 1.23462e+13         1      -1775.77871       1.20919e+13  0.60106  9.86645e+02      14
6449155011847092940      -1775.19009 1.23532e+13        -1      -1775.77871       1.21016e+13  0.60106  9.87239e+02      15
4593097647669372429      -1774.80408 7.33253e+12        -1      -1775.00607       2.64681e+12  0.98707  1.91944e-04       2
5143640288307395270      -1774.80552 7.08975e+12        -1      -1775.00790       2.54223e+12  0.98563  1.96258e-04       0
5090142258360873861      -1775.17161 1.46023e+13        -1      -1775.77871       1.41248e+13  0.61954  5.70910e+02      20
168126849921176544      -1774.80552 7.23585e+12        -1      -1775.00790       2.58683e+12  0.98563  2.00307e-04       2
1873538378006629199      -1774.80407 7.42065e+12        -1      -1775.00607       2.66689e+12  0.98707  1.94240e-04       1
8609431500298462695      -1774.80235 7.22446e+12        -1      -1775.00586       2.61357e+12  0.98880  1.76922e-04       0
```

Energies are in eV, prefactors and rates in s⁻¹.
The rate is the harmonic transition state theory rate, the prefactor times {math}`\exp(-E_\text{barrier}/kT)`, with {math}`kT` the thermal energy at 300 K.
`product` is the number of the product state, or -1 for a state not visited yet.

```{figure} ../fig/akmc/akmc_pt_processes.png
---
alt: Reactant, saddle and product energies of the 8 processes, with the process taken in red
align: center
---
Energy at the reactant, saddle and product of each process.
```

Three processes have barriers between 0.601 and 0.620 eV and end 0.012 eV above the reactant.
The other five have barriers between 0.985 and 0.989 eV and end about 0.78 eV above the reactant.
Their rates run from 1.76922e-04 to 2.00307e-04 s⁻¹, against 571 to 987 s⁻¹ for the low barriers.

## The KMC step

The step is recorded in `dynamics.txt`:

```{code-block} text
 step-number   reactant-id    process-id    product-id     step-time    total-time       barrier          rate        energy
-----------------------------------------------------------------------------------------------------------------------------
           0             0  1905950264122117683             1  1.329980e-03  1.329980e-03      0.601060  9.866450e+02  -1775.791149
```

The server took process 1905950264122117683 and advanced the clock by 1.33 ms.
Its reactant, saddle and product are in `states/0/procdata/`, named with the process number:

```{figure} ../fig/akmc/akmc_pt_process_taken.png
---
alt: Top view of the island before, at and after the saddle of the process taken
class: full-width
align: center
---
Top view of the process taken.
The island atoms are orange and the top surface layer is gray.
```

All seven island atoms move between 1.576 and 1.593 Å from reactant to product, and no other atom moves more than 0.144 Å.
The island slides as one piece.
Further server calls search for processes from state 1.

`docs/source/fig/akmc/akmc_pt.py` reruns this example and redraws the plots, and `akmc_cycle.dot` next to it is the source of the diagram.
