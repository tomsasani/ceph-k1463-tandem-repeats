import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

plt.rc("font", size=12)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]


ours = pd.read_csv(
    "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/chm13v2.0_maskedY_rCRS.palladium-v1.0.trgt.annotations.bed.gz",
    sep="\t",
    # nrows=1e5,
)
theirs = pd.read_csv(
    "data/gangstr.trgt_formatted.bed.gz",
    sep="\t",
    # nrows=1e5,
    names=["#chrom", "start", "end", "info"]
)
ours["who"] = "ours"
theirs["who"] = "theirs"

combined = pd.concat([ours, theirs])
combined["length"] = combined["end"] - combined["start"]

colors = sns.color_palette("colorblind", 2)

f, ax = plt.subplots()
bins = np.linspace(combined["length"].min(), combined["length"].max(), 1_000)
for g, gdf in combined.groupby("who"):
    lengths = gdf["length"].values
    total = np.sum(lengths)
    n50 = total / 2
    sorted_lengths = np.sort(lengths)[::-1]
    sorted_cumsum = np.cumsum(sorted_lengths)
    at_thresh = np.where(sorted_cumsum >= n50)[0][0]
    n50_len = sorted_lengths[at_thresh]
    hist, edges = np.histogram(lengths, bins=bins)
    hist_frac = hist / np.sum(hist)
    cumsum = np.cumsum(hist_frac)
    ax.plot(
        edges[:-1],
        cumsum,
        label=(
            f"this study (N50 = {n50_len}bp)"
            if g == "ours"
            else f"Mousavi et al. (N50 = {n50_len}bp)"
        ),
        c=colors[0] if g == "ours" else colors[1],
        zorder=0,
    )
    ax.axvline(x=n50_len, ls=":", c=colors[0] if g == "ours" else colors[1], zorder=-1)
ax.legend(shadow=True, title="Catalog")
ax.set_xscale("log")
ax.set_xlabel("TR locus size in reference genome (bp)")
ax.set_ylabel("Cumulative fraction of TR loci in catalog")
sns.despine(ax=ax)
f.savefig("plots/catalog_comparison.png", dpi=200)
