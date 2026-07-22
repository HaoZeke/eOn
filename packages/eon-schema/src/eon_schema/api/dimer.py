"""Dimer / min-mode API composition (L2)."""

from __future__ import annotations

from typing import Any, Optional, Union

from eon_schema._deps import ensure_import
from eon_schema.fields.enums import Accelerant, MinModeMethod

_DimerSpecCls: Any = None


def _dimer_spec_cls() -> Any:
    global _DimerSpecCls
    if _DimerSpecCls is not None:
        return _DimerSpecCls
    pd = ensure_import(
        "pydantic",
        install_spec="pydantic>=2",
        extra_hint="Or: pip install 'eon-schema[pydantic]'",
    )
    BaseModel = pd.BaseModel
    ConfigDict = pd.ConfigDict
    Field = pd.Field
    model_validator = pd.model_validator

    class DimerSpec(BaseModel):
        """Validated min-mode request for pyeonclient ``Dimer(...)``.

        Default method is improved. ``accelerant=gp`` only with improved.
        """

        model_config = ConfigDict(extra="forbid", use_enum_values=False)

        method: MinModeMethod = Field(default=MinModeMethod.improved)
        accelerant: Accelerant = Field(default=Accelerant.none)

        @model_validator(mode="after")
        def _gp_requires_improved(self) -> Any:
            if (
                self.accelerant is Accelerant.gp
                and self.method is not MinModeMethod.improved
            ):
                raise ValueError(
                    'accelerant="gp" is only valid with method="improved" '
                    f'(got method="{self.method.value}"). '
                    "GP-dimer is the improved-dimer + GP stack."
                )
            return self

        def core_kwargs(self) -> dict[str, str]:
            acc = (
                ""
                if self.accelerant is Accelerant.none
                else self.accelerant.value
            )
            return {"method": self.method.value, "accelerant": acc}

    _DimerSpecCls = DimerSpec
    return _DimerSpecCls


def __getattr__(name: str) -> Any:
    if name == "DimerSpec":
        return _dimer_spec_cls()
    raise AttributeError(name)


def parse_dimer_spec(
    method: Union[str, MinModeMethod] = "improved",
    accelerant: Union[str, Accelerant, None] = None,
    *,
    spec: Any = None,
) -> Any:
    if spec is not None:
        return spec
    acc: Any = Accelerant.none if accelerant in (None, "") else accelerant
    return _dimer_spec_cls().model_validate(
        {"method": method, "accelerant": acc}
    )
