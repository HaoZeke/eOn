"""Named force backends for Matter — simple key → factory map.

Register once, look up by name. No Matter C++ pollution; ASE / metatomic /
RGPOT factories live here as pure Python.

Usage::

    import pyeonclient as pyec
    from pyeonclient.backends import make_backend, list_backends

    pot = make_backend("ase", calculator=my_calc)
    pot = make_backend("metatomic", model_path="pet-mad.pt")
    pot = make_backend(
        "rgpot_metatomic",
        model_path="pet-mad.pt",
        engine_path="libmetatomic_engine.so",
    )
    m = pyec.Matter(pot, pyec.Parameters())
"""

from __future__ import annotations

from typing import Any, Callable

BackendFactory = Callable[..., Any]

# name → factory(**kwargs) → Potential
_REGISTRY: dict[str, BackendFactory] = {}


def register(name: str) -> Callable[[BackendFactory], BackendFactory]:
    """Decorator: ``@register("ase")`` adds a backend factory."""

    def deco(fn: BackendFactory) -> BackendFactory:
        key = name.strip().lower()
        if not key:
            raise ValueError("backend name must be non-empty")
        _REGISTRY[key] = fn
        return fn

    return deco


def list_backends() -> list[str]:
    """Sorted registered backend keys."""
    return sorted(_REGISTRY)


def make_backend(name: str, **kwargs: Any) -> Any:
    """Build a pyeonclient ``Potential`` by registry key."""
    key = name.strip().lower()
    try:
        fn = _REGISTRY[key]
    except KeyError as e:
        raise KeyError(
            f"unknown backend {name!r}; known: {list_backends()}"
        ) from e
    return fn(**kwargs)


# --- built-in factories -------------------------------------------------


@register("ase")
def _ase(*, calculator: Any, **_: Any) -> Any:
    """ASE Calculator wrapped as eOn Potential (no compile-time -Dwith_ase)."""
    if calculator is None:
        raise ValueError("backend 'ase' requires calculator=")
    from pyeonclient._core import make_potential_from_ase

    return make_potential_from_ase(calculator)


@register("metatomic")
@register("metatomic_fat")
def _metatomic_fat(
    *,
    model_path: str,
    device: str = "cpu",
    params: Any = None,
    **_: Any,
) -> Any:
    """Fat-linked ``PotType.METATOMIC`` (needs ``-Dwith_metatomic=true``)."""
    import pyeonclient as pyec

    if not model_path:
        raise ValueError("backend 'metatomic' requires model_path=")
    p = params if params is not None else pyec.Parameters()
    p.potential = pyec.PotType.METATOMIC
    p.metatomic_model_path = str(model_path)
    p.metatomic_device = str(device)
    return pyec.make_potential(p)


@register("rgpot_metatomic")
@register("metatomic_dlopen")
def _metatomic_dlopen(
    *,
    model_path: str,
    engine_path: str = "",
    device: str = "cpu",
    params: Any = None,
    **_: Any,
) -> Any:
    """RGPOT + ``backend=metatomic`` (dlopen ``libmetatomic_engine.so``)."""
    import pyeonclient as pyec

    if not model_path:
        raise ValueError("backend 'rgpot_metatomic' requires model_path=")
    p = params if params is not None else pyec.Parameters()
    p.potential = pyec.PotType.RGPOT
    p.rgpot_backend = "metatomic"
    p.rgpot_model_path = str(model_path)
    p.rgpot_device = str(device)
    if engine_path:
        p.rgpot_engine_path = str(engine_path)
    return pyec.make_potential(p)


@register("lj")
def _lj(*, params: Any = None, **_: Any) -> Any:
    """Built-in LJ (for smoke tests without models)."""
    import pyeonclient as pyec

    p = params if params is not None else pyec.Parameters()
    p.potential = pyec.PotType.LJ
    return pyec.make_potential(p)


__all__ = [
    "register",
    "list_backends",
    "make_backend",
]
