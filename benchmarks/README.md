# eOn ASV benchmarks

Python ASV suite that times the `eonclient` binary (saddle search, point
energy, …). See `bench_eonclient.py` and `data/`.

## Local

```bash
pixi install
# build + install eonclient into the pixi env (same as CI)
pixi run meson setup bbdir --prefix=$CONDA_PREFIX --libdir=lib --buildtype release
pixi run meson install -C bbdir
pixi run bash -c 'pip install asv asv-tachyon'
pixi run asv machine --yes
pixi run asv run -E "existing:$(which python)" --quick
pixi run asv publish
pixi run asv-tachyon install .asv/html
pixi run asv-tachyon serve .asv/html --open
```

## Dashboard CI

`.github/workflows/ci_asv_dashboard.yml` runs on `main`, weekly, and
`workflow_dispatch`. It publishes an [asv-tachyon](https://github.com/HaoZeke/asv_tachyon)
static site to the `gh-pages` branch under **`/bench/`**.

- Fork: `https://haozeke.github.io/eOn/bench/`
- Upstream docs site: `https://eondocs.org/bench/` (once merged)
- Optional custom host: `https://bench.eondocs.org` → CNAME to the Pages
  host that serves this tree (see workflow header comments)

PR comments still use [asv-perch](https://github.com/HaoZeke/asv-perch)
via `ci_benchmark.yml` + `ci_bench_commenter.yml`.
