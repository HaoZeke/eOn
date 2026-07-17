import pytest
from pydantic import ValidationError

from eon_schema.api import DimerSpec
from eon_schema.fields import Accelerant, MinModeMethod


def test_enums():
    assert MinModeMethod.improved.value == "improved"
    assert Accelerant.gp.value == "gp"


def test_dimer_spec():
    s = DimerSpec()
    assert s.method is MinModeMethod.improved
    s2 = DimerSpec(method=MinModeMethod.improved, accelerant=Accelerant.gp)
    assert s2.core_kwargs()["accelerant"] == "gp"
    with pytest.raises(ValidationError):
        DimerSpec(method="classic", accelerant="gp")
