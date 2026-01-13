import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.patches as patches


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


mutations = pd.read_csv(snakemake.input.mutations, sep="\t", dtype={"paternal_id": str, "sample_id": str})

# read in akshay decomposition
akshay = pd.read_csv(snakemake.input.akshay, sep="\t")
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


print (single_motif_mutations.shape, multi_motif_mutations.shape)


combined = pd.concat([single_motif_mutations, multi_motif_mutations])

combined["motif_change"] = combined["likely_denovo_size"] // combined["motif_size"]
combined["motif_change_leftover"] = combined["likely_denovo_size"] % combined["motif_size"]


f, axarr = plt.subplots(3, 2, figsize=(8, 8), sharey="row", sharex=True)
motif2i = dict(zip(range(1, 7), [0, 0, 1, 1, 2, 2]))
motif2j = dict(zip(range(1, 7), [0, 1, 0, 1, 0, 1]))

strs = combined[combined["motif_size"] <= 6]
strs = strs[strs["likely_denovo_size"].between(-20, 20)]
bins = np.arange(-20, 21, 1)

palette = sns.color_palette("colorblind", 2)

for m, mdf in strs.groupby("motif_size"):
    i, j = motif2i[m], motif2j[m]
    ax = axarr[i, j]
    hist, edges = np.histogram(mdf["likely_denovo_size"].values, bins=bins)
    exact_multiple_idxs = np.where(edges[:-1] % m == 0)[0]
    other_idxs = np.setdiff1d(np.arange(edges[:-1].shape[0]), exact_multiple_idxs)

    hist_frac = hist / np.sum(hist)

    rect1 = patches.Rectangle(
        (0, 0),
        max(edges),
        0.5,
        facecolor="green",
        zorder=-1,
        alpha=0.1,
    )
    rect2 = patches.Rectangle(
        (0, 0),
        min(edges),
        0.5,
        facecolor="red",
        zorder=-1,
        alpha=0.1,
    )

    ax.add_patch(rect1)
    ax.add_patch(rect2)


    if m == 1:
        for side in ("con", "exp"):
            arrow_start, arrow_tip = 5, 15
            if side == "con":
                arrow_start *= -1
                arrow_tip *= -1
            arr = patches.FancyArrowPatch(
                (arrow_start, 0.3),
                (arrow_tip, 0.3),
                arrowstyle="-|>",
                color="k",
                mutation_scale=20,
            )
            ax.add_patch(arr)
            lab = 5 if side == "exp" else -15
            ax.annotate("expansion" if side == "exp" else "contraction", xy=(lab, 0.2))

    for match, idxs in zip((True, False), (exact_multiple_idxs, other_idxs)):

        ind = edges[idxs]
        vals = hist_frac[idxs]
        lw = 1
        ec = "k" if match else "w"
        color = palette[1 - int(match)]
        ax.bar(
            ind,
            vals,
            1,
            lw=0.5,
            ec="w",
            color=color,
            label="perfect" if match else "imperfect",
            zorder=0,
        )
        title = f"Motif = {m} bp" if i == j == 0 else f"{m} bp"
        title += f"\n(n = {mdf.shape[0]})"
        ax.set_title(title)
        if i == 2:
            ax.set_xlabel("Inferred DNM size (bp)")
        if j == 0:
            ax.set_ylabel("Fraction of DNMs")
        if i == 2 and j == 0:
            ax.legend(shadow=True)
        for xpos in range(-20, 21):
            if xpos % m == 0 and m != 1:
                ax.axvline(xpos, ls="--", c="gainsboro", alpha=0.5, zorder=-1)
    sns.despine(ax=ax)
f.tight_layout()
f.savefig(snakemake.output.by_motif, dpi=200)


combined = combined[(combined["motif_change_leftover"] == 0) & (combined["motif_change"] != 0)]

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
    f, axarr = plt.subplots(combined["TR type"].nunique(), 1, figsize=(6, 8), sharex=True, sharey=False)

phase2idx = dict(zip(["dad", "mom", "unknown"], [0, 1, 2]))
type2idx = {"homopolymer": 0, "non-homopolymer STR": 1, "VNTR": 2}
colors = sns.color_palette("colorblind", 3)

val = "motif_change"

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
    sizes["frac"] = sizes["count"] / sizes["count"].sum()

    n_expansions = sub_df.query("motif_change > 0").shape[0]
    n_contractions = sub_df.query("motif_change < 0").shape[0]
    # res = ss.binomtest(n_expansions, n_expansions + n_contractions, alternative="less")
    # p = round(res.pvalue, 3)
    # print (sub_cols, p)

    ax = axarr[j]

    ind = sizes[val].values
    if SPLIT_BY_PHASE:
        ax.bar(ind, sizes["frac"].values, 1, ec="w", lw=1, color=colors[i])
        sns.despine(ax=axarr[i, j])
        ax.set_title(f"{tr_type}")
    else:
        ax.bar(ind, sizes["frac"].values, 1, ec="w", lw=1, color="darkgrey")
        rect1 = patches.Rectangle(
            (0, 0),
            20,
            sizes["frac"].max() * 1.1,
            facecolor="green",
            zorder=-1,
            alpha=0.1,
        )
        rect2 = patches.Rectangle(
            (0, 0),
            -20,
            sizes["frac"].max() * 1.1,
            facecolor="red",
            zorder=-1,
            alpha=0.1,
        )

        ax.add_patch(rect1)
        ax.add_patch(rect2)

        if tr_type == "homopolymer":
            for side in ("con", "exp"):
                arrow_start, arrow_tip = 5, 15
                if side == "con":
                    arrow_start *= -1
                    arrow_tip *= -1
                arr = patches.FancyArrowPatch(
                    (arrow_start, 0.125),
                    (arrow_tip, 0.125),
                    arrowstyle="-|>",
                    color="k",
                    mutation_scale=20,
                )
                ax.add_patch(arr)
                lab = 6 if side == "exp" else -14
                ax.annotate("expansion" if side == "exp" else "contraction", xy=(lab, 0.1))


        sns.despine(ax=ax)
        ax.set_title(f"{tr_type} (n = {sub_df.shape[0]})")

        if j == 1:
            ax.set_ylabel("Fraction of DNMs")
        if j == 2:
            ax.set_xlabel("Inferred DNM size (motif units)")

f.tight_layout()
f.savefig(snakemake.output.by_tr_type, dpi=200)
