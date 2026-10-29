import numpy as np
import pytest

gt = pytest.importorskip("genome_tools")

from ..api import load_data_results, save_data_results  # noqa: E402
from ..config import ModelConfig  # noqa: E402
from ..results import fit_footprints  # noqa: E402
from ..summaries import Summaries  # noqa: E402


def test_save_load_and_slice(tmp_path):
    from genome_tools.plotting.modular_plot.api import DataBundle

    rng = np.random.default_rng(0)
    x, v = rng.normal(-0.5, 0.6, (3, 40)), np.full((3, 40), 0.02)
    data = DataBundle(interval=gt.GenomicInterval("chr1", 1000, 1040))
    data.summaries = Summaries(("a", "b", "c"), x, v, v + 0.4, ("s1", "s2"), ("a", "b"), x[:2], v[:2])
    data.footprints = fit_footprints(x, v, v + 0.4, data.summaries.group_names, ModelConfig(n_draws=100))
    save_data_results(data, tmp_path / "dhs.npz")

    back = load_data_results(tmp_path / "dhs.npz")
    assert back.footprints.calls == data.footprints.calls
    assert back.summaries.group_names == ("a", "b", "c")

    sub = load_data_results(tmp_path / "dhs.npz", DataBundle(interval=gt.GenomicInterval("chr1", 1010, 1030)))
    assert sub.footprints.p_pattern.shape == (19,)
    assert np.allclose(sub.footprints.p_more, data.footprints.p_more[10:30])
    assert np.allclose(sub.summaries.x, x[:, 10:30])
