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
