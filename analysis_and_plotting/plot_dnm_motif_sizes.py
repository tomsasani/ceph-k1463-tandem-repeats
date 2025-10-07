import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import scipy.stats as ss
import tqdm


plt.rc("font", size=14)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

def get_size_in_motifs(row):
    size = row["likely_denovo_size"]
    motif_size = row["min_motiflen"]

    if size % motif_size != 0:
        return np.nan
    else:
        return size // motif_size

ASSEMBLY = "GRCh38"

mutations = pd.read_csv(f"{ASSEMBLY}.filtered.tsv", sep="\t", dtype={"paternal_id": str, "sample_id": str})

# read in akshay decomposition
akshay = pd.read_csv("data/K1463.CHM13v2.DNMs.416.demintr.output", sep="\t")
n_events = akshay.groupby("ID").size().reset_index().rename(columns={0: "n_events"})
akshay = akshay.merge(n_events)

akshay = akshay[akshay["n_events"] == 1]

akshay["akshay_motif"] = akshay["Info"].apply(lambda i: i.split(";")[5].split(":")[1] if len(i.split(";")[5].split(":")) > 1 else "UNK")
akshay["akshay_motif_size"] = akshay["Info"].apply(lambda i: len(i.split(";")[5].split(":")[1]) if len(i.split(";")[5].split(":")) > 1 else "UNK")
akshay = akshay[akshay["akshay_motif"] != "UNK"]

# get mutations with a single motif, these are easy
single_motif_mutations = mutations[mutations["n_motifs"] == 1]

# merge multi-motif mutations with akshays file to decompose
multi_motif_mutations = mutations[mutations["n_motifs"] > 1]
multi_motif_mutations = akshay.merge(multi_motif_mutations, left_on="ID", right_on="trid")

# subset to relevant columns
single_motif_mutations = single_motif_mutations[
    ["sample_id", "likely_denovo_size", "motifs", "min_motiflen", "phase"]
].rename(columns={"motifs": "motif", "min_motiflen": "motif_size"})

multi_motif_mutations = multi_motif_mutations[
    ["sample_id", "likely_denovo_size", "akshay_motif", "akshay_motif_size", "phase"]
].rename(columns={"akshay_motif": "motif", "akshay_motif_size": "motif_size"})

print (multi_motif_mutations)

combined = pd.concat([single_motif_mutations, multi_motif_mutations])
# combined = multi_motif_mutations

# combined = combined[~combined["sample_id"].str.startswith("200")]

combined["motif_change"] = combined["likely_denovo_size"] // combined["motif_size"]
combined["motif_change_leftover"] = combined["likely_denovo_size"] % combined["motif_size"]


f, axarr = plt.subplots(3, 2, figsize=(8, 8), sharey=True, sharex=True)
motif2i = dict(zip(range(1, 7), [0, 0, 1, 1, 2, 2]))
motif2j = dict(zip(range(1, 7), [0, 1, 0, 1, 0, 1]))

strs = combined[combined["motif_size"] <= 6]
strs = strs[strs["likely_denovo_size"].between(-20, 20)]
bins = np.arange(-20, 21, 1)

palette = sns.color_palette("colorblind", 2)

for m, mdf in strs.groupby("motif_size"):
    i, j = motif2i[m], motif2j[m]
    hist, edges = np.histogram(mdf["likely_denovo_size"].values, bins=bins)
    exact_multiple_idxs = np.where(edges[:-1] % m == 0)[0]
    other_idxs = np.setdiff1d(np.arange(edges[:-1].shape[0]), exact_multiple_idxs)

    hist_frac = hist / np.sum(hist)

    for match, idxs in zip((True, False), (exact_multiple_idxs, other_idxs)):

        ind = edges[idxs]
        vals = hist_frac[idxs]
        lw = 1
        ec = "k" if match else "w"
        color = palette[1 - int(match)]
        axarr[i, j].bar(ind, vals, 1, lw=0.5, ec="w", color=color, label="perfect" if match else "imperfect", zorder=0)
        title = f"Motif = {m} bp" if i == j == 0 else f"{m} bp"
        title += f"\n(n = {mdf.shape[0]})"
        axarr[i, j].set_title(title)
        if i == 2:
            axarr[i, j].set_xlabel("Inferred DNM size (bp)")
        if j == 0:
            axarr[i, j].set_ylabel("Fraction of DNMs")
        if i == 2 and j == 0:
            axarr[i, j].legend(shadow=True)
        for xpos in range(-20, 21):
            if xpos % m == 0 and m != 1:
                axarr[i, j].axvline(xpos, ls="--", c="gainsboro", alpha=0.5, zorder=-1)
    sns.despine(ax=axarr[i, j])
f.tight_layout()
f.savefig("strs.png", dpi=200)



combined = combined[(combined["motif_change_leftover"] == 0) & (combined["motif_change"] != 0)]
print (combined.shape)

# figure out "true" TR type
combined["TR type"] = combined["motif_size"].apply(
    lambda m: (
        "homopolymer" if m == 1 else "non-homopolymer STR" if 2 <= m <= 6 else "VNTR"
    )
)


SPLIT_BY_PHASE = False

if SPLIT_BY_PHASE:
    f, axarr = plt.subplots(combined["phase"].nunique(), combined["TR type"].nunique(), figsize=(16, 8), sharex=True, sharey="row")
else:
    f, axarr = plt.subplots(combined["TR type"].nunique(), 1, figsize=(6, 8), sharex=True)

phase2idx = dict(zip(["dad", "mom", "unknown"], [0, 1, 2]))
type2idx = dict(zip(combined["TR type"].unique(), range(combined["TR type"].nunique())))

colors = sns.color_palette("colorblind", 3)

val = "motif_change"
# val = "likely_denovo_size"

combined = combined[combined[val].between(-20, 20)]

group_cols = ["TR type", "phase"] if SPLIT_BY_PHASE else ["TR type"]

for (sub_cols), sub_df in combined.groupby(group_cols):
    if SPLIT_BY_PHASE:
        tr_type, phase = sub_cols
        i, j = phase2idx[phase], type2idx[tr_type]
    else:
        tr_type = sub_cols[0]
        phase = "all"
        i, j = 0, type2idx[tr_type]

    sizes = sub_df.groupby(val).size().reset_index().rename(columns={0: "count"})

    n_expansions = sub_df.query("motif_change > 0").shape[0]
    n_contractions = sub_df.query("motif_change < 0").shape[0]
    res = ss.binomtest(n_expansions, n_expansions + n_contractions, alternative="less")
    p = round(res.pvalue, 3)
    print (sub_cols, p)

    ind = sizes[val].values
    if SPLIT_BY_PHASE:
        axarr[i, j].bar(ind, sizes["count"].values, 1, ec="w", lw=1, color=colors[i])
        sns.despine(ax=axarr[i, j])
        axarr[i, j].set_title(f"{tr_type}")
    else:
        axarr[j].bar(ind, sizes["count"].values, 1, ec="w", lw=1, color="darkgrey")
        sns.despine(ax=axarr[j])
        axarr[j].set_title(f"{tr_type} (n = {sub_df.shape[0]})")

        if j == 1:
            axarr[j].set_ylabel("Number of DNMs")
        if j == 2:
            axarr[j].set_xlabel("Inferred size of DNM (motif units)")

f.tight_layout()
f.savefig("sizes.motifs.png", dpi=200)
