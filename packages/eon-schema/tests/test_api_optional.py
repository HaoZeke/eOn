import pytest

def test_enums_no_pydantic():
    from eon_schema.fields import MinModeMethod, Accelerant
    assert MinModeMethod.improved.value == "improved"
    assert Accelerant.gp.value == "gp"

def test_dimer_spec_needs_pydantic():
    pytest.importorskip("pydantic")
    from eon_schema.api import DimerSpec
    from eon_schema.fields import Accelerant, MinModeMethod
    from pydantic import ValidationError

    s = DimerSpec()
    assert s.method is MinModeMethod.improved
    s2 = DimerSpec(method=MinModeMethod.improved, accelerant=Accelerant.gp)
    assert s2.core_kwargs()["accelerant"] == "gp"
    with pytest.raises(ValidationError):
        DimerSpec(method="classic", accelerant="gp")
