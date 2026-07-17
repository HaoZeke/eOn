import pytest
from pydantic import ValidationError

from eon_schema.api import DimerSpec, NebSpec
from eon_schema.fields import Accelerant, MinModeMethod, PathInit


def test_enums():
    assert MinModeMethod.improved.value == "improved"
    assert Accelerant.gp.value == "gp"
    assert PathInit.file.value == "file"


def test_dimer_spec():
    s = DimerSpec()
    assert s.method is MinModeMethod.improved
    s2 = DimerSpec(method=MinModeMethod.improved, accelerant=Accelerant.gp)
    assert s2.core_kwargs()["accelerant"] == "gp"
    with pytest.raises(ValidationError):
        DimerSpec(method="classic", accelerant="gp")


def test_neb_spec_cookbook_defaults():
    s = NebSpec(
        n_images=10,
        path_init=PathInit.file,
        path_list="idppPath.dat",
        energy_weighted=True,
        ci_mmf=True,
        max_iterations=1000,
        force_tolerance=0.01,
        random_seed=706253457,
    )
    assert s.n_images == 10
    assert s.path_list == "idppPath.dat"
    assert s.energy_weighted is True
    assert s.ci_mmf is True
    assert s.core_accelerant() == ""
