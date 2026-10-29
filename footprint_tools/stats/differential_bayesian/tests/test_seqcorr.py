import numpy as np

from footprint_tools.stats.differential_bayesian import seqcorr


def test_recovers_planted_sequence_effect():
    rng = np.random.default_rng(0)
    stats, n = [], 300
    for k in range(12):
        seq = "".join(rng.choice(list("ACGT"), n + 2 * seqcorr.W + 1))
        effect = np.array([1.0 if seq[seqcorr.W + j:seqcorr.W + j + 2] == "CG" else 0.0 for j in range(n)])
        exp = np.full(n, 200.0)
        obs = rng.poisson(exp * 2 ** (effect - effect.mean())).astype(float)
        stats.append((f"chr{k % 3}", obs, exp, seq))
    coef, r2 = seqcorr.fit(stats)
    assert r2 > 0.5
    x = stats[0][3]
    cg = np.array([x[seqcorr.W + j:seqcorr.W + j + 2] == "CG" for j in range(n)])
    corr = seqcorr.correction(x, coef)
    assert corr[cg].mean() - corr[~cg].mean() > 0.7
