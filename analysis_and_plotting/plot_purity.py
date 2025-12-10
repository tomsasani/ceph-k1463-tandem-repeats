import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as ss
import matplotlib.patches as patches

plt.rc("font", size=11)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

rng = np.random.default_rng(42)

mutations = pd.read_csv(snakemake.input.mutations, sep="\t", dtype={"sample_id": str})

for c in mutations.columns:
   if np.any(mutations[c].isna()): print (c)

mutations = mutations[mutations["haplotype_in_parent_consensus"] != "UNK"]

mutations["precursor_AL"] = mutations["precursor_sequence_in_parent"].apply(lambda s: len(s))
mutations["untransmitted_AL"] = mutations["untransmitted_sequence_in_parent"].apply(lambda s: len(s))


mutations["direction"] = mutations["likely_denovo_size"].apply(
    lambda d: "expansion" if d > 0 else "contraction" if d < 0 else "interruption"
)


f, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4), sharey=True)

palette = sns.color_palette("colorblind", 2)

vals = []

for i, (val, ax) in enumerate(zip(("AL", "AP"), (ax1, ax2))):

    tidy = mutations.melt(
        id_vars=["sample_id", "trid", "genotype", "likely_denovo_size"],
        value_vars=[f"precursor_{val}", f"untransmitted_{val}"],
        var_name="allele_status",
        value_name="size",
    )

    tidy["is_precursor"] = tidy["allele_status"].apply(lambda s: 1 if "precursor" in s else 0)

    # wilcoxon signed rank test on allele purities. alternative hypothesis
    # of "left" indicates that the precursor APs are on the whoel lower than the untransmitted
    stat, p = ss.wilcoxon(
        mutations[f"precursor_{val}"].values,
        y=mutations[f"untransmitted_{val}"].values,
        alternative="greater",
    )

    a, b = mutations[f"precursor_{val}"].values, mutations[f"untransmitted_{val}"].values
    diffs = a - b

    bin_edges = np.linspace(min(diffs), max(diffs), 50)

    v_hist, v_edges = np.histogram(diffs, bins=bin_edges)
    v_hist_frac = v_hist / np.sum(v_hist)
    v_hist_cumsum = np.cumsum(v_hist_frac)

    ind = np.arange(v_edges[1:].shape[0])

    ax.bar(v_edges[1:], v_hist, v_edges[1:] - v_edges[:-1], ec="w", color="k", lw=1,  zorder=1)
    rect1 = patches.Rectangle((0, 0), max(v_edges), 1e3, facecolor='green', zorder=-1, alpha=0.15)
    rect2 = patches.Rectangle((0, 0), min(v_edges), 1e3, facecolor='red', zorder=-1, alpha=0.15)

    # annotate with arrows
    for side in ("less", "more"):
        
        arrow_start = max(v_edges) / 4 if side == "more" else min(v_edges) / 4
        arrow_tip = arrow_start * 3 
        arr = patches.FancyArrowPatch((arrow_start, 1e2), (arrow_tip, 1e2),
                               arrowstyle='-|>', color="k", mutation_scale=20)
        ax.add_patch(arr)
        lab = (arrow_tip + arrow_start) / 3. if side == "more" else (arrow_tip + arrow_start) / 1.3

        annotation = "transmitted allele\n    is"
        if val == "AP":
            annotation += f" {side} pure"
        else:
            annotation += " longer" if side == "more" else " shorter"
       

    # Add the patch to the axes
    ax.add_patch(rect1)
    ax.add_patch(rect2)
    ax.set_title(f"Wilcoxon p = {p:.2e}")
    ax.set_xlabel("Difference between allele lengths..." if val == "AL" else "...and allele purities")
    sns.despine(ax=ax)
    ax.set_yscale("log")

ax1.set_ylabel(f"Number of TR DNMs")

f.tight_layout()
f.savefig(snakemake.output.png, dpi=200)
