"""Process id allocation must not freeze on duplicate processtable rows.

get_num_procs() is len(distinct ids). After a duplicate id appears, that
length freezes and reusing it as the next id overwrites procdata. New
registrations must use get_next_process_id() == max(id)+1 (eOn-cj3o).
"""

from __future__ import annotations

import os
from unittest import mock

import pytest

from eon.akmcstate import AKMCState


def _make_akmc_state(tmp_path, table_body: str = ""):
    from eon import fileio as io

    state_dir = tmp_path / "state_0"
    state_dir.mkdir()
    (state_dir / "procdata").mkdir()
    proctable = state_dir / "processtable"
    proctable.write_text(AKMCState.processtable_header + table_body)

    state = object.__new__(AKMCState)
    state.config = mock.Mock()
    state.statelist = mock.Mock()
    state.statelist.kT = 0.025
    state.path = str(state_dir)
    state.number = 0
    state.procs = None
    state.proc_repeat_count = None
    state.procdata_path = os.path.join(state.path, "procdata")
    state.reactant_path = os.path.join(state.path, "reactant.con")
    state.proctable_path = str(proctable)
    state.search_result_path = os.path.join(state.path, "search_results.txt")
    state.tar_path = os.path.join(state.path, "procdata.tar")
    state.info = io.ini(os.path.join(state.path, "info"))
    state.info.set("MetaData", "kT", 0.025)
    return state


# Shape of the Cu V/SIA observation: 8 rows, 3 distinct ids (0,1,2), counter
# freezes at 3 if next id is taken from len(procs).
_DUP_TABLE = """\
      0          1.00000 1.00000e+12         1          0.50000       1.00000e+12  0.74164  3.47609e-01       0
      1          1.10000 1.00000e+12         2          0.60000       1.00000e+12  0.70427  1.47483e+00       0
      2          1.20000 1.00000e+12         3          0.70000       1.00000e+12  0.70086  1.68270e+00       0
      2          0.90000 1.00000e+12         4          0.10000       1.00000e+12  0.10270  1.88217e+10       0
      2          1.15000 1.00000e+12         5          0.65000       1.00000e+12  0.68879  2.68439e+00       0
      2          0.91000 1.00000e+12         6          0.11000       1.00000e+12  0.10171  1.95599e+10       0
      2          1.14000 1.00000e+12         7          0.64000       1.00000e+12  0.65673  9.27689e+00       0
      2          0.92000 1.00000e+12         8          0.12000       1.00000e+12  0.10271  1.88176e+10       0
"""


def test_get_num_procs_is_distinct_count_not_row_count(tmp_path):
    state = _make_akmc_state(tmp_path, _DUP_TABLE)
    assert state.get_num_procs() == 3
    assert sorted(state.procs.keys()) == [0, 1, 2]
    # Last row for id 2 wins.
    assert state.procs[2]["barrier"] == pytest.approx(0.10271)


def test_get_next_process_id_is_max_plus_one_not_len(tmp_path):
    state = _make_akmc_state(tmp_path, _DUP_TABLE)
    # Old bug: get_num_procs() == 3 would reuse id 2 and clobber procdata.
    assert state.get_num_procs() == 3
    assert state.get_next_process_id() == 3


def test_next_id_after_empty_table_is_zero(tmp_path):
    state = _make_akmc_state(tmp_path, "")
    assert state.get_num_procs() == 0
    assert state.get_next_process_id() == 0


def test_append_refuses_duplicate_id(tmp_path):
    state = _make_akmc_state(tmp_path, _DUP_TABLE)
    state.load_process_table()
    with pytest.raises(RuntimeError, match="refusing to clobber process id 2"):
        state.append_process_table(
            id=2,
            saddle_energy=1.0,
            prefactor=1e12,
            product=9,
            product_energy=0.0,
            product_prefactor=1e12,
            barrier=0.5,
            rate=1.0,
            repeats=0,
        )


def test_append_with_next_id_does_not_clobber(tmp_path):
    state = _make_akmc_state(tmp_path, _DUP_TABLE)
    nid = state.get_next_process_id()
    assert nid == 3
    state.append_process_table(
        id=nid,
        saddle_energy=1.0,
        prefactor=1e12,
        product=99,
        product_energy=0.0,
        product_prefactor=1e12,
        barrier=0.2,
        rate=1e5,
        repeats=0,
    )
    assert state.get_num_procs() == 4
    assert state.get_next_process_id() == 4
    assert 2 in state.procs
    assert state.procs[3]["product"] == 99
    text = open(state.proctable_path).read()
    assert any(line.split()[:1] == ["3"] for line in text.splitlines()[1:])


def test_force_reload_sees_external_rewrite(tmp_path):
    row0 = (
        "      0          1.00000 1.00000e+12         1          0.50000"
        "       1.00000e+12  0.50000  1.00000e+00       0\n"
    )
    row5 = (
        "      5          1.00000 1.00000e+12         2          0.50000"
        "       1.00000e+12  0.40000  1.00000e+00       0\n"
    )
    state = _make_akmc_state(tmp_path, row0)
    assert state.get_next_process_id() == 1
    # External rewrite (e.g. amsel KDB replay) with a higher max id.
    open(state.proctable_path, "w").write(
        AKMCState.processtable_header + row0 + row5
    )
    assert state.get_next_process_id() == 6
