from pathlib import Path


ROOT = Path(__file__).parents[1]
LOCK_FILES = sorted((ROOT / "condaEnvs").glob("*.conda-lock.yml"))


def test_pixi_lock_records_plotting_stack_floor():
    text = (ROOT / "pixi.lock").read_text(encoding="utf-8")

    assert "rgpycrumbs>=1.10.4,<2" in text
    assert "rgpycrumbs>=1.3" not in text


def test_conda_locks_record_plotting_and_readcon_floors():
    assert LOCK_FILES

    for path in LOCK_FILES:
        text = path.read_text(encoding="utf-8")
        assert "rgpycrumbs: '>=1.10.4,<2" in text, path
        assert "rgpycrumbs: '>=1.3" not in text, path
        assert "readcon-chemfiles: ==0.14.7" in text, path
        assert "readcon-chemfiles: ==0.14.5" not in text, path


def test_rgpot_wrap_uses_the_abi_compatible_revision():
    text = (ROOT / "subprojects/rgpot.wrap").read_text(encoding="utf-8")

    assert "revision = b19227f9197a4c927633f1225f1e7546794a6f1e" in text
    assert "revision = 473033de21a26dc25adb73df3b3640eed0dc0074" not in text
