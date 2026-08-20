"""Monorepo: JobResult Cap'n Proto file + results.dat adapters."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_schema_file_in_monorepo():
    p = ROOT / "schema" / "eon_job_result.capnp"
    assert p.is_file()
    text = p.read_text(encoding="utf-8")
    assert "struct JobResult" in text
    assert "struct JobRequest" in text
    assert "struct Geometry" in text
    assert "struct OptimizerProvenance" in text
    assert "struct EngineCompatibility" in text
    assert "struct LandfoldArtifact" in text
    assert "landfoldArtifacts @31 :List(LandfoldArtifact);" in text
    # Flat geometry, not con text blobs
    assert "positions @0" in text
    assert ".con" not in text or "not .con text" in text


def test_adapters_importable_from_package():
    import sys

    sys.path.insert(0, str(ROOT / "packages" / "eon-schema" / "src"))
    from eon_schema.jobs import job_result_capnp_path, results_dat_to_dict

    assert job_result_capnp_path().is_file()
    d = results_dat_to_dict("0 termination_reason\ngood termination_reason_text\n")
    assert d["termination_reason"] == 0


def test_results_dat_adapter_exposes_optimizer_provenance():
    import sys

    sys.path.insert(0, str(ROOT / "packages" / "eon-schema" / "src"))
    from eon_schema.jobs import job_result_scalars_from_results_dat

    parsed = job_result_scalars_from_results_dat(
        """0 termination_reason
minimization job_type
xtsci optimizer_backend
eon.optimizer.v1 optimizer_provenance_schema
eon.compatibility.v1 compatibility_schema
eon engine_id
2.11.1 engine_version
abc123 engine_build_identity
eon.objective compatibility_engine_protocol_family
1 compatibility_engine_protocol_major
0 compatibility_engine_protocol_minor
1 compatibility_engine_abi_major
0 compatibility_engine_abi_minor
2 compatibility_engine_layout_revision
1 optimizer_xts_abi_major
0 optimizer_xts_abi_minor
2 optimizer_xts_abi_layout
"""
    )
    assert parsed["optimizer"]["backend"] == "xtsci"
    assert parsed["optimizer"]["schema"] == "eon.optimizer.v1"
    assert parsed["optimizer"]["xts_abi"] == {"major": 1, "minor": 0, "layout": 2}
    assert parsed["compatibility"]["engine"] == {
        "id": "eon",
        "version": "2.11.1",
        "build_identity": "abc123",
    }
    assert parsed["compatibility"]["engine_compatibility"] == {
        "schema": "eon.compatibility.v1",
        "engineId": "eon",
        "protocolFamily": "eon.objective",
        "protocolMajor": 1,
        "protocolMinor": 0,
        "abiMajor": 1,
        "abiMinor": 0,
        "layoutRevision": 2,
        "buildIdentity": "abc123",
    }


def test_results_dat_adapter_exposes_compatibility_provenance():
    import sys

    sys.path.insert(0, str(ROOT / "packages" / "eon-schema" / "src"))
    from eon_schema.jobs import job_result_scalars_from_results_dat

    parsed = job_result_scalars_from_results_dat(
        """0 termination_reason
lj potential_type
eon.compatibility.v1 compatibility_schema
3 compatibility_readcon_spec_version
0.14.7 compatibility_readcon_min_version
2.11.1 engine_version
abc123 engine_build_identity
"""
    )
    assert parsed["compatibility"]["schema"] == "eon.compatibility.v1"
    assert parsed["compatibility"]["readcon"] == {
        "spec_version": 3,
        "min_version": "0.14.7",
    }
    assert parsed["compatibility"]["engine"] == {
        "id": "lj",
        "version": "2.11.1",
        "build_identity": "abc123",
    }


def test_results_dat_adapter_does_not_upgrade_unknown_compatibility_schema():
    import sys

    sys.path.insert(0, str(ROOT / "packages" / "eon-schema" / "src"))
    from eon_schema.jobs import job_result_scalars_from_results_dat

    parsed = job_result_scalars_from_results_dat(
        """eon.compatibility.v0 compatibility_schema
eon engine_id
eon.objective compatibility_engine_protocol_family
1 compatibility_engine_protocol_major
0 compatibility_engine_protocol_minor
1 compatibility_engine_abi_major
0 compatibility_engine_abi_minor
2 compatibility_engine_layout_revision
abc123 engine_build_identity
"""
    )
    assert "engine_compatibility" not in parsed["compatibility"]
