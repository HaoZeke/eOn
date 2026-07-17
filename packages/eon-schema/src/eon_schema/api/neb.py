"""NEB API composition (L2)."""

from __future__ import annotations

from typing import Any, Union

from eon_schema._deps import ensure_import
from eon_schema.fields.enums import Accelerant, PathInit

_NebSpecCls: Any = None


def _neb_spec_cls() -> Any:
    global _NebSpecCls
    if _NebSpecCls is not None:
        return _NebSpecCls
    pd = ensure_import(
        "pydantic",
        install_spec="pydantic>=2",
        extra_hint="Or: pip install 'eon-schema[pydantic]'",
    )
    BaseModel = pd.BaseModel
    ConfigDict = pd.ConfigDict
    Field = pd.Field

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
            try:
                from pyeonclient._core import NEBInit  # type: ignore
            except ImportError as e:
                raise ImportError(
                    "NebSpec.apply_to_parameters requires pyeonclient"
                ) from e
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

    _NebSpecCls = NebSpec
    return _NebSpecCls


def __getattr__(name: str) -> Any:
    if name == "NebSpec":
        return _neb_spec_cls()
    raise AttributeError(name)


def parse_neb_spec(
    accelerant: Union[str, Accelerant, None] = None,
    *,
    path_init: Union[str, PathInit] = "idpp",
    n_images: int = 10,
    spec: Any = None,
) -> Any:
    if spec is not None:
        return spec
    acc: Any = Accelerant.none if accelerant in (None, "") else accelerant
    return _neb_spec_cls().model_validate(
        {"accelerant": acc, "path_init": path_init, "n_images": n_images}
    )
