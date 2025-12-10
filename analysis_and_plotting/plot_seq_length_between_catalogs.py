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

f, ax = plt.subplots()
bins = np.linspace(combined["length"].min(), combined["length"].max(), 1_000)
for g, gdf in combined.groupby("who"):
    lengths = gdf["length"].values
    hist, edges = np.histogram(lengths, bins=bins)
    hist_frac = hist / np.sum(hist)
    cumsum = np.cumsum(hist_frac)
    ax.plot(edges[1:], cumsum, label="this study" if g == "ours" else "GangSTR - Mousavi et al. (2019)")
ax.legend(shadow=True, title="Catalog")
ax.set_xscale("log")
ax.set_xlabel("TR locus size in reference genome (bp)")
ax.set_ylabel("Fraction of TR loci in catalog")
sns.despine(ax=ax)
f.savefig("o.png", dpi=200)
