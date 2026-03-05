import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.patches as mpatches

plt.rc("font", size=11, )

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

TRID = "chr8_2623352_2623487_trsolve"

PALETTE = "colorblind"

ceph_motifs = pd.read_csv(
    f"trviz/ceph/{TRID}.GRCh38.key.tsv",
    sep="\t",
    names=["sequence", "encoded_motif", "n_obs"],
)["sequence"].unique()


key = pd.read_csv(f"trviz/hprc/{TRID}.GRCh38.key.tsv", sep="\t", names=["sequence", "encoded_motif", "n_obs"])
seq = pd.read_csv(f"trviz/hprc/{TRID}.GRCh38.encoded.tsv", sep="\t")

seq["allele_length"] = seq["encoded_motifs"].apply(lambda s: len(s))
seq["allele_length"] = seq["original_sequence"].apply(lambda s: len(s))
# seq["contains_motif"] = seq["encoded_motifs"].str.contains(MOTIF)
seq["sample"] = seq["sample_id"].apply(lambda s: s.split("_")[0])
seq["hap"] = seq["sample_id"].apply(lambda s: int(s.split("_")[6]))


# for each "encoded motif" in the key...
res = []
for encoded_motif in key["encoded_motif"].unique():
    # figure out how many copies of that encoded motif every sample had
    motif_counts_per_sample = []
    for i, row in seq.iterrows():
        sample_id = row["sample_id"]
        sample_encoded_motifs = row["encoded_motifs"]
        motif_count = sum([m == encoded_motif for m in list(sample_encoded_motifs)])
        res.append(
            {
                "sample_id": sample_id,
                "encoded_motif": encoded_motif,
                "motif_count": motif_count,
            }
        )
res_df = pd.DataFrame(res)

# get the motifs present in the CEPH individuals
# key = pd.read_csv(f"trviz/hprc/{TRID}.{GENOME}.key.tsv", sep="\t", names=["motif", "encoded_motif", "count"])
encoded2motif = dict(zip(key["encoded_motif"], key["sequence"]))

ceph_motifs = [
    "GAGGCGCCAGGAGAGCGCT",
    "GAGGCGCCAGGAGAGAGCGCT",
    "GAGCGCCAAGCGCT",
    "GGAGGCGCCAGGAGAGAGCGCT",
    "GGAGGCGCCAGGAGAGCGCT",
    "GAGGCGCCAGGAGCGCGCT",
    "GGAGGCGCCAGGAGCGCGCT",
]

res_df["motif"] = res_df["encoded_motif"].apply(lambda m: encoded2motif[m])

res_df = res_df[res_df["motif"].isin(ceph_motifs)]

# identify the unique characters in the encoded
# motif arrays, including "-" gaps
uniq_chars = res_df["motif"].unique()
uniq_vals = np.arange(len(uniq_chars))
cmap = sns.color_palette(PALETTE, res_df["motif"].nunique())

f, (ax, ax1) = plt.subplots(1, 2, figsize=(3, 4))
bins = np.arange(1, res_df["motif_count"].max() + 1)
for m, mdf in res_df.groupby('encoded_motif'):
    vals = mdf["motif_count"].values
    hist, edges = np.histogram(vals, bins=bins)
    ax.hist(vals, histtype="step")
    ax1.hist(vals, histtype="step")

f.savefig("hprc_hist.png")

f, (ax, ax1) = plt.subplots(2, 1, figsize=(5, 3), height_ratios=[1, 3])

for _ax in (ax, ax1):
    sns.stripplot(
        data=res_df,
        x="encoded_motif",
        y="motif_count",
        hue="motif",
        # marker_kws={"edgecolor": "w", "linewidth": 1},
        palette=PALETTE,
        edgecolor="w",
        linewidth=0.5,
        ax=_ax,
        legend=False,
    )
ax.set_ylabel("")
ax1.set_xlabel("Motif")

ax.set_ylim(73, 77)  # outliers only
ax1.set_ylim(-2, 28)  # most of the data

# hide the spines between ax and ax2
ax.tick_params(bottom=False, labelbottom=False)
sns.despine(ax=ax, bottom=True)
sns.despine(ax=ax1)

d = .015  # how big to make the diagonal lines in axes coordinates
# # arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d*3, +d*3), **kwargs)        # top-left diagonal

kwargs.update(transform=ax1.transAxes)  # switch to the bottom axes
ax1.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal

# make custom legend
# get colors used in palette
custom_legend = []
for char, color in zip(uniq_chars, cmap):
    custom_legend.append(mpatches.Patch(color=color, lw=1, ec="w"))


ax1.set_ylabel("Motif count on HPRC alleles")
ax.set_xlabel("")
f.tight_layout()
f.savefig("hprc_motifs.png", dpi=200)