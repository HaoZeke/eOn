"""results.dat timing footer must be value-key for parse_results."""
from __future__ import annotations

from eon.fileio import parse_results
from eon_schema.jobs import job_result_scalars_from_results_dat, results_dat_to_dict


def test_parse_results_value_key_timing():
    text = (
        "0 termination_reason\n"
        "good termination_reason_text\n"
        "1.234567 time_seconds\n"
        "0.5 user_time\n"
        "0.1 system_time\n"
    )
    d = parse_results(__import__("io").StringIO(text))
    assert abs(d["time_seconds"] - 1.234567) < 1e-12
    assert abs(d["user_time"] - 0.5) < 1e-12
    assert abs(d["system_time"] - 0.1) < 1e-12


def test_inverted_key_value_timing_does_not_yield_float_keys():
    """Historical inverted footer put the key first; parsers must not
    treat that as the value-key contract."""
    inverted = "time_seconds 1.234567\n"
    d = parse_results(__import__("io").StringIO(inverted))
    # key becomes the number string; no float under time_seconds
    assert "time_seconds" not in d or not isinstance(d.get("time_seconds"), float)


def test_job_result_adapter_reads_timing():
    text = (
        "0 termination_reason\n"
        "2.5 time_seconds\n"
        "1.0 user_time\n"
        "0.25 system_time\n"
    )
    mapped = job_result_scalars_from_results_dat(text)
    assert abs(mapped["wall_time_seconds"] - 2.5) < 1e-12
    assert abs(mapped["user_time_seconds"] - 1.0) < 1e-12
    assert abs(mapped["system_time_seconds"] - 0.25) < 1e-12
    assert abs(results_dat_to_dict(text)["time_seconds"] - 2.5) < 1e-12
