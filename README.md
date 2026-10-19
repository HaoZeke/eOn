# asv-results (orphan branch)

Long-lived [airspeed velocity](https://github.com/airspeed-velocity/asv) result
JSON for eOn. Updated by the **ASV dashboard** GitHub Actions workflow.

## Layout

```
.asv/results/
  <machine>/
    *.json
  benchmarks.json   # when present
```

## Rules

- **Orphan branch** — no shared history with `main`.
- Only benchmark outputs (and this README). No source trees.
- Force-pushed as a whole tree by CI after each successful run (safe: single writer).

## Consumers

- CI restores this tree before `asv run` so history accumulates.
- `asv publish` + `asv-tachyon install` builds the static dashboard under `gh-pages/bench/`.
