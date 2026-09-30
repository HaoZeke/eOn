---
myst:
  html_meta:
    "description": "Tunnelling splittings between two minima in eOn, and the thermal rate through a saddle below the crossover: WKB along a NEB band and the ring-polymer instanton job."
    "keywords": "eOn instanton, tunnelling splitting, two-level system, WKB, ring polymer, instanton rate."
---

# Tunnelling splittings and rates

A two-level system (TLS) in a glass is a pair of adjacent minima that the
structure tunnels between at about one kelvin. Its tunnelling splitting
{math}`\Delta_0` and asymmetry {math}`\Delta` set the TLS energy
{math}`E = \sqrt{\Delta^2 + \Delta_0^2}`. eOn estimates {math}`\Delta_0` in two
ways:

| | NEB band (WKB) | `job = instanton` |
|---|---|---|
| Path | the minimum energy path | the path of least imaginary-time action |
| Dimensions | one, along the mass-weighted band | every free degree of freedom |
| Cost | the NEB itself | one batch of forces per iteration over the beads, plus a Hessian per bead |
| Output | first frame of `neb.con` | `instanton.con` and `results.dat` |

With `mode = rate` the same job estimates the thermal rate through a saddle
instead of a splitting. That path is below.

Energies are in eV, lengths in Å and masses in amu throughout, so
mass-weighted lengths are in amu^0.5 Å.

## WKB along a NEB band

Every nudged elastic band (NEB) job writes `reaction_coordinate_mw`, the mass-weighted arc length, on
each frame of `neb.con`. The first frame also carries the band's Wentzel-Kramers-Brillouin (WKB) estimate.
The keys are:

| Key | Meaning |
|---|---|
| `hbar_omega_reactant`, `hbar_omega_product` | Wells: {math}`\hbar\omega` of each well along the band, from a fit of {math}`a s^2 + b s^3` to the images within half the barrier |
| `tunnel_action` | Action: {math}`S = \hbar^{-1} \int \sqrt{2 (V(s) - E)}\, ds` over the forbidden region |
| `tunnel_splitting` | Estimate: {math}`\Delta_0 = (\hbar\omega / \pi) e^{-S}`, with {math}`\omega` the geometric mean of the wells |
| `tls_energy` | Energy: {math}`\sqrt{\Delta^2 + \Delta_0^2}` |
| `tunnel_deep_wells` | Flag: 1 when both barriers exceed {math}`\hbar\omega`; below that, WKB is the wrong tool |

The profile between images is a monotone cubic, so it cannot dip below the
data. The level {math}`E` is the higher of the two harmonic ground states. A
structure without masses, or a band whose end is flat, leaves these keys out.

WKB along the band is exact in one dimension up to its semiclassical error.
When the path curves, the tunnelling cuts the corner, and the transverse
zero-point energy changes along the way. In the two-dimensional test valley
below, both effects together put the band estimate a factor of 2.8 below the
exact splitting.

## The instanton job

The instanton is a path of {math}`P + 1` beads in imaginary time
{math}`\beta\hbar`, with its ends fixed at the two minima. It minimises the
discretised Euclidean action

```{math}
S = \sum_j \frac{|q_{j+1} - q_j|^2}{2\,\delta\tau} + \delta\tau \sum_j V(q_j),
\qquad \delta\tau = \beta\hbar / P,
```

in mass-weighted coordinates {math}`q`. The splitting comes from the ratio of
the off-diagonal to the diagonal imaginary-time propagator, both taken in the
same steepest-descent approximation.

```{math}
\Delta_0 = 2\hbar \sqrt{\frac{S_0}{2\pi\hbar\,\delta\tau}}
\sqrt{\frac{\det J_\text{well}}{\det' J}}\; e^{-(S - S_\text{well})/\hbar}.
```

Here {math}`J` is the Hessian of {math}`S` over the interior beads. The
prime leaves out its zero mode, the kink's position in imaginary time.
{math}`S_0 = \int |\dot q|^2 d\tau`, and {math}`J_\text{well}` is the same
Hessian with every bead at a minimum.

```{code-block} ini
[Main]
job = instanton

[Potential]
potential = rgpot

[Instanton]
reactant_filename = reactant.con
product_filename = product.con
; start from a converged band instead of the straight line
initial_path = neb.con
beads = 256
beta_hbar_omega = 30
force_tolerance = 1e-3
; one finite-difference Hessian every 4 beads, linear in between
hessian_stride = 4
```

`beta_hbar_omega` sets the imaginary time in units of {math}`1/\omega` of
the stiffer minimum along the line between them. It must be long enough for
the kink to relax into both wells. `results.dat` reports
`instanton_mode_separation`, the second eigenvalue of {math}`J` over its zero
mode. Values above {math}`10^3` mean the kink is isolated; small values mean
{math}`\beta\hbar` is too short. With no atom fixed, the product is aligned
to the reactant first: its mass-weighted mean displacement is removed, and
for a cluster its best rotation as well.

With that reactant rotation removed, each iteration evaluates every interior bead of the kink in one call. Under
`[RgpotPot] ranks_per_image`, that call spreads the beads over the Car-Parrinello molecular dynamics (CPMD)
calculator groups the same way a NEB spreads its images.

`instanton.con` writes one frame per bead, with `imaginary_time_fs`. The keys are:

| Key | Meaning |
|---|---|
| `tunnel_splitting_instanton` | Splitting: {math}`\Delta_0`, eV |
| `instanton_action` | Action: {math}`(S - S_\text{well})/\hbar` |
| `tls_energy_instanton` | Energy: {math}`\sqrt{\Delta^2 + \Delta_0^2}`, eV |
| `tunnel_asymmetry` | Asymmetry: {math}`V(\text{product}) - V(\text{reactant})`, eV |
| `instanton_temperature_K` | Temperature: {math}`1/(k_B \beta)` for the imaginary time used |
| `instanton_mode_separation` | Separation: how well the kink's translation separates from the other modes |
| `instanton_symmetric` | Symmetry: 1 when {math}`\beta|\Delta| < 0.1` |
| `instanton_beta_asymmetry` | Magnitude: {math}`\beta|\Delta|` |

The propagator ratio measures the splitting {math}`\Delta_0` when the two wells lie within
a small fraction of {math}`k_B T` of each other. `instanton_symmetric = 0`
flags a pair outside that window. The job still writes the path and the
action, but no `tunnel_splitting_instanton`, and reports success: the flag
says why. For such a pair, set `mode = rate`, give the saddle and a
temperature below the crossover, and read the rate in the section below.

## Which path object

NEB images and ring-polymer beads are different objects. An image is a point
on a path in configuration space between two minima. Its springs are
fictitious. The parallel force is removed. A bead is
one imaginary-time slice of a single quantum system. Its springs are
physical, with stiffness fixed by the temperature and the number of beads.

The columns are:

| Object | Points | Springs | What it returns |
|---|---|---|---|
| `mode = splitting` | open string between two minima, started from a band when one is present | Euclidean action, no tangent projection | tunnelling splitting when the wells are close in energy |
| `mode = rate` | closed ring through one saddle | same action, stiffness set by {math}`T` and {math}`N` | thermal rate below the crossover temperature |
| Centroid potential of mean force (PMF) | one ring per image, centroid held on the image | sampled, not minimised | quantum free-energy barrier along the path |
| Harmonic centroid string | a ring at each image, optimised | local harmonic quantum correction | a free-energy estimate as good as that harmonic well |

The first two rows are `job = instanton`. The centroid potential of mean
force is a constrained path-integral molecular dynamics sample, one
thermostatted ring per image. That sample is not an optimisation, and it
does not belong in this job. A string of harmonically corrected rings is a
different calculation again, and this page does not implement it.

## The instanton rate

`mode = rate` finds the ring-polymer instanton for the thermal rate out of
the reactant through a first-order saddle, at a temperature below the
crossover {math}`T_c = \hbar\omega_b / (2\pi k_B)`, with {math}`\omega_b`
the imaginary frequency at the saddle (Richardson and Althorpe 2009;
Richardson 2016). The instanton is a first-order saddle of the discretised
Euclidean action on a closed path of {math}`N` beads, the formulation of
Einarsdóttir et al. (2012) and Ásgeirsson, Arnaldsson and Jónsson (2018):

```{math}
U_N = \sum_j V(q_j) + \sum_j \frac{|q_{j+1} - q_j|^2}{2\beta_N^2\hbar^2},
\qquad \beta_N = \beta / N, \qquad q_N = q_0,
```

and the rate is

```{math}
k\,Z_r = \frac{1}{\beta_N\hbar}\sqrt{\frac{B_N}{2\pi\beta_N\hbar^2}}
\;\prod_k{}' \frac{1}{\beta_N\hbar|\omega_k|}\; e^{-\beta_N U_N},
```

with {math}`B_N = \sum_j |q_{j+1} - q_j|^2`, {math}`\omega_k^2` the
eigenvalues of the ring Hessian of {math}`U_N` without its zero mode (the
ring's translation in imaginary time), and {math}`Z_r` the ring-polymer
partition function of the harmonic reactant. Rigid-body modes drop out of
both.

```{code-block} ini
[Main]
job = instanton

[Instanton]
mode = rate
reactant_filename = reactant.con
saddle_filename = saddle.con
; a band over the barrier seeds the ring and gives the WKB rate along it
initial_path = neb.con
beads = 64
temperatures = 300, 250, 200, 150
hessian_final = recompute
hessian_stride = 4
```

The search is Newton eigenvector following on the ring Hessian, block
cyclic tridiagonal in the beads, with the bead blocks seeded from the
saddle's Hessian and Bofill-updated from the gradient differences; each
step costs one batch of {math}`N/2 + 1` force calls, the ring held
symmetric under imaginary-time reversal. On the one-dimensional Eckart
barrier this converges in 4 to 7 steps from either seed. The determinant,
its inertia and the linear solves go through a block LU of the open chain
plus a low-rank correction for the closure and the zero mode, so the
{math}`Nf \times Nf` matrix is never formed; the cost is
{math}`O(N f^3)` in {math}`f` degrees of freedom.

The ring starts, in this order of preference, from the ring of the previous
temperature (`temperatures` runs from the highest down, each ring seeding
the next), from the path in `initial_path` mapped onto imaginary time by the
period condition {math}`\oint ds / \sqrt{2(V(s) - E)} = \beta\hbar`, or from
the saddle's unstable mode. `bead_ladder = true` converges a quarter of the
beads first, then half, then all, when no path seeds the ring.

The prefactor needs a Hessian at every bead. `hessian_final = recompute`
takes finite-difference Hessians on every `hessian_stride`-th bead of the
half ring and interpolates linearly in between, at `beads / (2
hessian_stride)` Hessians; `updated` keeps the Bofill-updated blocks the
search ends with, at no cost and lower accuracy.

Output, per temperature, is one row of `rate_instanton.dat` (T, {math}`T_c`,
beads, convergence, {math}`U_N`, negative modes, {math}`\ln(k / \mathrm{s}^{-1})`,
{math}`k`, harmonic TST, the effective barrier {math}`-k_B T \ln(2\pi\hbar\beta k)`,
and {math}`\ln k` from the one-dimensional WKB integral along the path when a
path was given), a frame per bead in `instanton.con` (the last temperature)
and `instanton_<T>K.con` when several temperatures ran. `results.dat`
carries the last temperature:

| Key | Meaning |
|---|---|
| `rate_instanton`, `rate_instanton_log` | {math}`k` in 1/s and {math}`\ln(k\,\mathrm{s})`; the rate itself underflows a double in deep tunnelling |
| `rate_htst`, `rate_htst_log` | classical harmonic transition-state theory at the same T |
| `rate_wkb_path_log` | the Kemble WKB rate along `initial_path` relative to the harmonic reactant well |
| `barrier_effective_instanton` | {math}`-k_B T \ln(2\pi\hbar\beta k)`, eV |
| `instanton_crossover_K`, `instanton_temperature_K` | {math}`T_c` and T |
| `instanton_negative_modes`, `instanton_zero_mode` | one, and a number near zero, for a converged ring |
| `instanton_ring_potential`, `instanton_bN` | {math}`U_N` in eV and {math}`B_N` in amu Å² |

A temperature at or above {math}`T_c` is reported and skipped: the ring
collapses onto the saddle and classical transition-state theory with a
quantum prefactor applies. A ring whose Hessian has a second negative mode
is written but carries no rate.

## Checks

These checks use two Catch2 cases. The curved-valley case is `Instanton splitting in a curved valley matches the exact gap`. The corner case is `The instanton cuts the corner the minimum energy path takes`.
The potential is {math}`V = V_0 (x^2 - 1)^2 + \tfrac{K}{2} (y - C (1 - x^2))^2` at unit mass
with {math}`K = 4` eV/Å². The exact gap comes from a fourth-order
finite-difference Hamiltonian, converged to {math}`10^{-5}`.

| {math}`V_0` / eV | {math}`C` / Å | exact {math}`\Delta_0` / eV | instanton / exact | WKB on the valley floor / exact |
|---|---|---|---|---|
| 0.12 | 0.35 | 2.520e-5 | 1.17 | 0.35 |
| 0.30 | 0.35 | 9.818e-8 | 1.08 | |
| 0.12 | 0 | 2.040e-5 | 1.12 | 1.02 |
| 0.30 | 0 | 1.195e-7 | 1.07 | |

The instanton's error falls as the barrier deepens, the regime glass TLS sit
in. The same cases tie the C++ path and splitting to an independent
implementation of the discretisation to {math}`2 \times 10^{-3}`.

For the rate, `The Eckart rate instanton matches the exact flux to its
semiclassical error` compares {math}`k Z_r` through the symmetric Eckart
barrier {math}`V_0 / \cosh^2(x/a)` ({math}`V_0 = 0.425` eV, {math}`a =
0.734` amu^0.5 Å, {math}`T_c = 150` K) with the exact flux
{math}`(2\pi\hbar)^{-1}\int P(E) e^{-\beta E} dE` from Eckart's transmission
probability, at {math}`T = 0.5\,T_c` and {math}`0.35\,T_c`:

| beads | instanton / exact |
|---|---|
| 64 | 0.94 to 0.96 |
| 128 | 0.93 to 0.94 |
| {math}`N \to \infty` (1/N² extrapolation) | 0.928 |

In one dimension the instanton is the steepest-descent evaluation of the
WKB thermal integral, so its limit shares the uniform WKB error; the Kemble
integral along the path gives the same 0.928. `The rate instanton of a cubic
well matches its decay rate` checks the metastable cubic well at
{math}`\beta\hbar\omega_0 = 30` against the Caldeira and Leggett
zero-temperature decay rate: 128 beads within 20 percent, 256 within 5, the
extrapolation within 1. `The ring spectrum from the block chain matches the
dense Hessian` ties the chain's determinant and inertia to a dense
eigendecomposition, and `A rigid mode leaves the instanton rate unchanged`
checks the rigid-mode bookkeeping.

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: INST_
keyprefix: inst-
---
```
