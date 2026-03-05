import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as ss
from FAILING_TRIDS import FAIL_RECURRENTS

plt.rc("font", size=14)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]


rng = np.random.default_rng(42)


mutations = pd.read_csv("csv/filtered_for_plots/CHM13v2.TOPUP.tsv", sep="\t")

recurrents = pd.read_csv("csv/recurrent/CHM13v2.TOPUP.recurrent.tsv", sep="\t")
recurrents = recurrents[recurrents["sufficient_cohort_depth"] == 1]
recurrents = recurrents[~recurrents["trid"].isin(FAIL_RECURRENTS)]



genome = pd.read_csv("data/sequence_report.tsv", sep="\t")

acros = [3, 14, 15, 21, 22]
acros = [f"chr{c}" for c in acros]

mutations["acro"] = mutations["#chrom"].isin(acros)
# mutations = mutations[mutations["acro"] == False]

# mutations = mutations[mutations["phase"] != "unknown"]

chrom2size = dict(zip(genome["UCSC style name"], genome["Seq length"]))

mutations["mid"] = (mutations["end"] + mutations["start"]) / 2.
mutations["rel_pos"] = mutations.apply(lambda row: row["mid"] / chrom2size[row["#chrom"]], axis=1)
palette = sns.color_palette("colorblind", 3)
bins = np.arange(0, 1, 0.01)
f, axarr = plt.subplots(3, 2, figsize=(12, 8), sharey="row")
n_bs = 10_000
for i, (m, mdf) in enumerate(mutations.groupby("simple_motif_size")):

    for j, (a, adf) in enumerate(mdf.groupby("acro")):
        ax = axarr[i, j]
        ax.axline((0, 0), slope=1, color='darkgrey', ls=":", zorder=-1)

        alab = " (acrocentric)" if a else ""
        vals = adf["rel_pos"].values
        bs = np.zeros((n_bs, bins.shape[0] - 1))
        for bi in range(n_bs):
            _vals = rng.choice(vals, replace=True, size=vals.shape[0])
            _hist, _ = np.histogram(_vals, bins=bins)
            _fracs = _hist / np.sum(_hist)
            _cumu = np.cumsum(_fracs)
            bs[bi] = _cumu

        hist, edges = np.histogram(vals, bins=bins)
        fracs = hist / np.sum(hist)
        cumu = np.cumsum(fracs)
        ax.plot(bins[:-1], cumu, color=palette[i], zorder=1)
        ci_u = np.percentile(bs, 97.5, axis=0)
        ci_l = np.percentile(bs, 2.5, axis=0)
        print (ci_u)
        print (ci_l)
        ax.fill_between(bins[:-1], ci_l, ci_u, color=palette[i], alpha=0.2, zorder=2)
        res = ss.ks_1samp(vals, ss.uniform.cdf, alternative="two-sided")
        p = round(res.pvalue, 4)

    
        sns.despine(ax=ax)
        ax.set_title(m + alab)
        if i == 2:
            ax.set_xlabel("Relative position on chromosome")
        if j == 0:
            ax.set_ylabel("Cumulative fraction\nof TR DNMs (+/- 95% CI)")
f.tight_layout()
f.savefig("hist.png", dpi=200)


# recurrents["mid"] = (recurrents["end"] + recurrents["start"]) / 2.
# recurrents["rel_pos"] = recurrents.apply(lambda row: row["mid"] / chrom2size[row["#chrom"]], axis=1)
# palette = sns.color_palette("colorblind", 3)
# bins = np.arange(0, 1, 0.01)
# f, ax = plt.subplots(figsize=(8, 6))
# ax.hist(recurrents["rel_pos"].values, bins=bins)

# ax.set_xlabel("Relative position on chromosome")
# ax.set_xlabel("Relative position on chromosome")
# ax.set_ylabel("Number of TR DNMs")
# f.tight_layout()
# f.savefig("hist.png", dpi=200)
