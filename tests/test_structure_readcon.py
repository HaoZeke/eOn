"""Structure ↔ readcon.ConFrame is the canonical I/O path (no ASE)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import readcon

from eon.structure import Structure, Atoms
from eon import fileio as io

DATA = Path(__file__).resolve().parent / "data"


def test_atoms_alias_is_structure():
    assert Atoms is Structure


def test_from_to_conframe_roundtrip_free_fixed():
    frame = readcon.ConFrame(
        cell=[10.0, 10.0, 10.0],
        angles=[90.0, 90.0, 90.0],
        atoms=[
            readcon.Atom("Cu", 0.0, 0.0, 0.0, [True, True, True], 1, 63.5),
            readcon.Atom("Cu", 1.5, 0.0, 0.0, [False, False, False], 2, 63.5),
        ],
        prebox_header=["t", ""],
    )
    s = Structure.from_conframe(frame)
    assert len(s) == 2
    assert s.free[0] == 0.0
    assert s.free[1] == 1.0
    assert s.names == ["Cu", "Cu"]
    back = s.to_conframe()
    assert list(back.atoms[0].fixed) == [True, True, True]
    assert list(back.atoms[1].fixed) == [False, False, False]
    np.testing.assert_allclose(
        [[a.x, a.y, a.z] for a in back.atoms], s.r, atol=1e-12
    )


def test_loadcon_returns_structure_via_readcon():
    path = DATA / "server" / "Pt_Heptamer_oneLayer" / "pos.con"
    s = io.loadcon(str(path))
    assert isinstance(s, Structure)
    ref = readcon.read_con(str(path))[0]
    assert len(s) == len(ref)
    assert s.names[0] == ref.atoms[0].symbol
    np.testing.assert_allclose(s.r[0], [ref.atoms[0].x, ref.atoms[0].y, ref.atoms[0].z])


def test_fileio_does_not_import_ase():
    import eon.fileio as fio
    import inspect

    src = inspect.getsource(fio)
    assert "import ase" not in src
    assert "from ase" not in src
    assert "import readcon" in src
