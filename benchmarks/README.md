# eOn ASV benchmarks

Python ASV suite that times the `eonclient` binary (saddle search, point
energy, …). See `bench_eonclient.py` and `data/`.

## Local

```bash
pixi install
pixi run meson setup bbdir --prefix=$CONDA_PREFIX --libdir=lib --buildtype release
pixi run meson install -C bbdir
pixi run bash -c 'pip install asv asv-tachyon'
pixi run asv machine --yes
pixi run asv run -E "existing:$(which python)" --quick
pixi run asv publish
pixi run asv-tachyon install .asv/html
pixi run asv-tachyon serve .asv/html --open
```

## History storage

Long-lived results live on the **orphan** git branch `asv-results`
(`.asv/results/**` only). The **ASV dashboard** workflow restores that branch
before each run and force-pushes updates after publish so history is not lost
to Actions cache eviction.

## Dashboard CI

`.github/workflows/ci_asv_dashboard.yml` runs on `main`, weekly, and
`workflow_dispatch`. It:

1. Restores `asv-results`
2. Builds `eonclient`, runs ASV, publishes with [asv-tachyon](https://github.com/HaoZeke/asv_tachyon)
3. Pushes results back to `asv-results`
4. Deploys UI to `gh-pages` → **`/bench/`**
5. Optionally deploys to Netlify (`NETLIFY_AUTH_TOKEN` + `NETLIFY_SITE_ID`) for **bench.eondocs.org**

| Surface | URL |
|---------|-----|
| Docs | https://eondocs.org/ |
| Bench (path) | https://eondocs.org/bench/ |
| Bench (host) | https://bench.eondocs.org/ |
| Fork preview | https://haozeke.github.io/eOn/bench/ |

PR comments still use [asv-perch](https://github.com/HaoZeke/asv-perch)
via `ci_benchmark.yml` + `ci_bench_commenter.yml`.
