"""The experimental namespace collects the non-stable APIs."""
import nsevt
import nsevt.experimental as ex


def test_experimental_exports_expected_symbols():
    assert set(ex.__all__) == {
        "block_conformal", "split_conformal", "ConformalBand",
        "twoscale_trend", "wasserstein_decomposition", "TwoScaleResult",
    }


def test_top_level_aliases_are_the_same_objects():
    # backward compatibility: the top-level names resolve to the same objects
    for name in ex.__all__:
        assert getattr(nsevt, name) is getattr(ex, name)
