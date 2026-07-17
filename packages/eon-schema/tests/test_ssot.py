from eon_schema.ssot import capnp_path, catalog_path, load_catalog, covered_sections

def test_capnp_exists():
    assert capnp_path().is_file()
    assert "MainOptions" in capnp_path().read_text(encoding="utf-8")

def test_catalog_loads():
    cat = load_catalog()
    assert "sections" in cat
    secs = covered_sections()
    assert len(secs) >= 1
