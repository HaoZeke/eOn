#!/usr/bin/env bash
# Compare setup notes: SocketNWChem vs RgpotPot(NWChem). Times are wall-clock
# for a single force call harness when an engine library is available.
set -euo pipefail
OUT="${1:-rgpot_nwchem_vs_socket.txt}"
{
  echo "=== Rgpot NWChem vs SocketNWChem (setup comparison) ==="
  echo "Date: $(date -Is)"
  echo
  echo "SocketNWChem: eOn listens; NWChem connects (i-PI). See"
  echo "  client/potentials/SocketNWChem/ and unit_tests/run_nwchem_test.sh"
  echo "  Requires nwchem on PATH and a .nwi with driver socket."
  echo
  echo "RgpotPot NWChem: in-process NWChemPot (dlopen libnwchemc), no potserv."
  echo "  eOn: potential=RGPOT, [RgpotPot] backend=NWChem -> capnp configure+calculate in-process"
  echo "  Engine: NWCHEMC_LIBRARY=... (or RGPOT_NWCHEMC_ENGINE / [RgpotPot] engine_library)"
  echo "  Fake CI path: rgpot's libnwchemc_fake_engine.so"
  echo
  if [[ -n "${NWCHEMC_LIBRARY:-}${RGPOT_NWCHEMC_ENGINE:-}" ]]; then
    echo "Engine library set (run test_rgpot_pot for timings)"
  else
    echo "Set NWCHEMC_LIBRARY and run: meson test -C <build> test_rgpot_pot --suite rgpot"
  fi
} | tee "$OUT"
