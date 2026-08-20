from eon_schema.config.models import OptimizerConfig, RefineConfig


def test_xtsci_is_an_explicit_optimizer_schema_value():
    assert OptimizerConfig(opt_method="xtsci").opt_method == "xtsci"
    assert RefineConfig(refine_opt_method="xtsci").refine_opt_method == "xtsci"
