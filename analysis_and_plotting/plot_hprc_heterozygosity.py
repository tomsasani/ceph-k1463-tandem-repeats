import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

plt.rc("font", size=12)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

df = pd.read_csv(
    "csv/hprc/combined.GRCh38.TOPUP.heterozygosity.tsv",
    sep="\t",
    names=[
        "chrom",
        "start",
        "end",
        "min_motiflen",
        "sample_id",
        "is_het",
        "is_denovo",
    ],
)
df = df[df["chrom"] != "chrX"]
df["motif"] = df["min_motiflen"].apply(lambda m: m if m < 7 else 7)
# get frac polymorphic at each TRID

df["is_denovo"] = df["is_denovo"].apply(lambda d: "yes" if d else "no")
# get number of HPRC samples with genotypes at each site
n_hprc = df.groupby(["chrom", "start", "end"]).agg(n_hprc=("sample_id", lambda s: len(set(s)))).reset_index()
df = df.merge(n_hprc)
# require all HPRC samples to be genotyped at each site
df = df.query("n_hprc >= 100")

grouped = df.groupby(["chrom", "start", "end", "motif", "is_denovo"]).agg(n=("is_het", "sum")).reset_index()
grouped["n"] = grouped["n"] / 100


f, ax = plt.subplots(
    figsize=(9, 6),
)
sns.violinplot(
    data=grouped,
    x="motif",
    y="n",
    hue="is_denovo",
    dodge=True,
    ax=ax,
    palette="colorblind",
    density_norm="width",
    inner=None,
    legend=True,
)
xlim = ax.get_xlim()
ylim = ax.get_ylim()
for vi, violin in enumerate(ax.collections):
    bbox = violin.get_paths()[0].get_extents()
    x0, y0, width, height = bbox.bounds
    # if it's the first of the hue, crop on the right.
    # if it's the second, crop on the left.
    mid_x = width / 2
    if vi % 2 == 1:
        violin.set_clip_path(
            plt.Rectangle((x0 + mid_x, y0), width / 2, height, transform=ax.transData),
        )
    else:
        violin.set_clip_path(
            plt.Rectangle((x0, y0), width / 2, height, transform=ax.transData),
        )

old_len_collections = len(ax.collections)

sns.stripplot(
    data=grouped,
    x="motif",
    y="n",
    hue="is_denovo",
    dodge=True,
    ax=ax,
    alpha=0.1,
    palette="colorblind",
    legend=False,
)
for di, dots in enumerate(ax.collections[old_len_collections:]):
    offset = -0.1 if di % 2 == 1 else 0.1
    dots.set_offsets(dots.get_offsets() + np.array([offset, 0]))
sns.despine(ax=ax)

ax.legend(
    title="Locus with DNM in K1463?",
    shadow=True,
    loc="upper center",
    bbox_to_anchor=(0.5, 1.15),
)
f.tight_layout()
ax.set_xlabel("Minimum motif size in locus")
ax.set_ylabel("Fraction of HPRC genotypes that are heterozygous")
ax.set_xticks(range(7))
ax.set_xticklabels(list(map(str, range(1, 7))) + ["7+"])
ax.set_yticks(np.arange(-0.2, 1.4, 0.2))
ax.set_yticklabels([""] + [str(round(x, 2)) for x in np.arange(0, 1.2, 0.2)] + [""])
f.savefig("hprc.png", dpi=200)
