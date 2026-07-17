"""Named force backends for Matter — simple key → factory map.

Register once, look up by name. No Matter C++ pollution; ASE / metatomic /
RGPOT factories live here as pure Python.

Usage::

    import pyeonclient as pyec
    from pyeonclient.backends import make_backend, list_backends

    pot = make_backend("ase", calculator=my_calc)
    pot = make_backend("metatomic", model_path="pet-mad.pt")
    pot = make_backend("ase_metatomic", model_path="pet-mad.pt")
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
_METATOMIC_LOAD_COMPAT_INSTALLED = False


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


def load_atomistic_model_compat(path: str, extensions_directory: Any = None) -> Any:
    """Load an exported metatomic model for ASE / Python use.

    Upstream ``load_atomistic_model`` re-wraps a scripted ``AtomisticModel`` and
    crashes on PET-MAD (and peers) when TorchScript does not expose
    ``_model_capabilities_outputs_names``.  Inject that list from
    ``capabilities().outputs`` before re-wrapping.
    """
    import torch
    from metatomic.torch import (
        AtomisticModel,
        check_atomistic_model,
        load_model_extensions,
    )

    path = str(path)
    load_model_extensions(
        path,
        str(extensions_directory) if extensions_directory is not None else None,
    )
    check_atomistic_model(path)
    model = torch.jit.load(path)
    try:
        model._model_capabilities_outputs_names  # noqa: B018
    except AttributeError:
        setattr(
            model,
            "_model_capabilities_outputs_names",
            list(model.capabilities().outputs.keys()),
        )
    return AtomisticModel(model, model.metadata(), model.capabilities())


def ensure_metatomic_load_compat() -> None:
    """Patch ``metatomic.torch.load_atomistic_model`` (and ASE calculator import).

    Safe to call repeatedly. Required so ``MetatomicCalculator(path)`` works for
    PET-MAD without the RecursiveScriptModule attribute error.
    """
    global _METATOMIC_LOAD_COMPAT_INSTALLED
    if _METATOMIC_LOAD_COMPAT_INSTALLED:
        return

    import metatomic.torch as mt
    import metatomic.torch.model as model_mod

    model_mod.load_atomistic_model = load_atomistic_model_compat
    mt.load_atomistic_model = load_atomistic_model_compat
    try:
        import metatomic_ase._calculator as calc_mod

        calc_mod.load_atomistic_model = load_atomistic_model_compat
    except ImportError:
        pass
    _METATOMIC_LOAD_COMPAT_INSTALLED = True


def make_metatomic_ase_calculator(
    model_path: str,
    *,
    device: str = "cpu",
    **calc_kwargs: Any,
) -> Any:
    """``metatomic_ase.MetatomicCalculator`` with PET-MAD-safe model load."""
    ensure_metatomic_load_compat()
    from metatomic_ase import MetatomicCalculator

    return MetatomicCalculator(str(model_path), device=device, **calc_kwargs)


# --- built-in factories -------------------------------------------------


@register("ase")
def _ase(*, calculator: Any, **_: Any) -> Any:
    """ASE Calculator wrapped as eOn Potential (no compile-time -Dwith_ase)."""
    if calculator is None:
        raise ValueError("backend 'ase' requires calculator=")
    from pyeonclient._core import make_potential_from_ase

    return make_potential_from_ase(calculator)


@register("ase_metatomic")
def _ase_metatomic(
    *,
    model_path: str = "",
    device: str = "cpu",
    calculator: Any = None,
    **calc_kwargs: Any,
) -> Any:
    """ASE ``MetatomicCalculator`` → eOn Potential (works with PET-MAD .pt)."""
    if calculator is None:
        if not model_path:
            raise ValueError(
                "backend 'ase_metatomic' requires model_path= or calculator="
            )
        calculator = make_metatomic_ase_calculator(
            model_path, device=device, **calc_kwargs
        )
    return _ase(calculator=calculator)


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
    "load_atomistic_model_compat",
    "ensure_metatomic_load_compat",
    "make_metatomic_ase_calculator",
]
