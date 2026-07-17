"""Pydantic models for the chemist-facing pyeonclient surface.

Prefer the shared package when installed::

    pip install 'eon-schema[pydantic]'

Otherwise fall back to in-tree models (same API). Neither path is required
for core Dimer/NEB: C++ + light validation always enforce
``accelerant=\"gp\"`` only with ``method=\"improved\"``.
"""

from __future__ import annotations

try:
    from eon_schema.api import DimerSpec, NebSpec, parse_dimer_spec, parse_neb_spec
    from eon_schema.fields import Accelerant, MinModeMethod, PathInit
except ImportError:
    # In-tree fallback (optional pydantic)
    from enum import Enum
    from typing import Any, Optional, Union

    try:
        from pydantic import BaseModel, ConfigDict, Field, model_validator
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "pyeonclient.models requires pydantic v2 or eon-schema[pydantic].\n"
            "  pip install 'pyeonclient[models]'\n"
            "  # or: pip install 'eon-schema[pydantic]'\n"
            "  # uv: uv pip install 'eon-schema[pydantic]'\n"
            "Core Dimer/NEB work without pydantic; C++ enforces "
            'accelerant=\"gp\" only with method=\"improved\".'
        ) from e

    class MinModeMethod(str, Enum):
        improved = "improved"
        classic = "classic"
        lanczos = "lanczos"
        davidson = "davidson"

    class Accelerant(str, Enum):
        none = "none"
        gp = "gp"

    class PathInit(str, Enum):
        linear = "linear"
        idpp = "idpp"
        idpp_collective = "idpp_collective"
        sidpp = "sidpp"
        sidpp_zbl = "sidpp_zbl"

    class DimerSpec(BaseModel):
        model_config = ConfigDict(extra="forbid", use_enum_values=False)
        method: MinModeMethod = Field(default=MinModeMethod.improved)
        accelerant: Accelerant = Field(default=Accelerant.none)

        @model_validator(mode="after")
        def _gp_requires_improved(self) -> "DimerSpec":
            if (
                self.accelerant is Accelerant.gp
                and self.method is not MinModeMethod.improved
            ):
                raise ValueError(
                    'accelerant="gp" is only valid with method="improved" '
                    f'(got method="{self.method.value}").'
                )
            return self

        def core_kwargs(self) -> dict[str, str]:
            acc = (
                ""
                if self.accelerant is Accelerant.none
                else self.accelerant.value
            )
            return {"method": self.method.value, "accelerant": acc}

    class NebSpec(BaseModel):
        model_config = ConfigDict(extra="forbid", use_enum_values=False)
        accelerant: Accelerant = Field(default=Accelerant.none)
        path_init: PathInit = Field(default=PathInit.idpp)
        n_images: int = Field(default=10, ge=1)

        def core_accelerant(self) -> str:
            return (
                ""
                if self.accelerant is Accelerant.none
                else self.accelerant.value
            )

        def apply_to_parameters(self, params: Any) -> Any:
            from pyeonclient._core import NEBInit

            mapping = {
                PathInit.linear: NEBInit.LINEAR,
                PathInit.idpp: NEBInit.IDPP,
                PathInit.idpp_collective: NEBInit.IDPP_COLLECTIVE,
                PathInit.sidpp: NEBInit.SIDPP,
                PathInit.sidpp_zbl: NEBInit.SIDPP_ZBL,
            }
            params.neb_init_method = mapping[self.path_init]
            params.neb_images = int(self.n_images)
            return params

    def parse_dimer_spec(
        method: str | MinModeMethod = "improved",
        accelerant: str | Accelerant | None = None,
        *,
        spec: Optional[DimerSpec] = None,
    ) -> DimerSpec:
        if spec is not None:
            return spec
        acc: str | Accelerant = (
            Accelerant.none if accelerant in (None, "") else accelerant
        )
        return DimerSpec.model_validate({"method": method, "accelerant": acc})

    def parse_neb_spec(
        accelerant: str | Accelerant | None = None,
        *,
        path_init: str | PathInit = "idpp",
        n_images: int = 10,
        spec: Optional[NebSpec] = None,
    ) -> NebSpec:
        if spec is not None:
            return spec
        acc: str | Accelerant = (
            Accelerant.none if accelerant in (None, "") else accelerant
        )
        return NebSpec.model_validate(
            {"accelerant": acc, "path_init": path_init, "n_images": n_images}
        )


__all__ = [
    "Accelerant",
    "DimerSpec",
    "MinModeMethod",
    "NebSpec",
    "PathInit",
    "parse_dimer_spec",
    "parse_neb_spec",
]
