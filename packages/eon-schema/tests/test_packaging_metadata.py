from pathlib import Path
import tomllib


def test_rgpycrumbs_extra_matches_plotting_stack_floor():
    package_root = Path(__file__).parents[1]
    metadata = tomllib.loads((package_root / "pyproject.toml").read_text())

    assert metadata["project"]["optional-dependencies"]["rgpycrumbs"] == [
        "rgpycrumbs>=1.10.4,<2"
    ]
