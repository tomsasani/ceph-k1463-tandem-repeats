import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import numpy as np

from matplotlib.colors import ListedColormap

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]
plt.rc("font", size=12)


mutations = pd.read_csv("GRCh38.filtered.tsv", sep="\t", dtype={"sample_id": str, "paternal_id": str})
print (mutations.shape)
# read in akshay decomposition
akshay = pd.read_csv("data/K1463.CHM13v2.DNMs.416.demintr.output", sep="\t", dtype={"sample_id": str})
n_events = akshay.groupby("ID").size().reset_index().rename(columns={0: "n_events"})
akshay = akshay.merge(n_events)

akshay = akshay[akshay["n_events"] == 1]

akshay["akshay_motif"] = akshay["Info"].apply(
    lambda i: (
        i.split(";")[5].split(":")[1] if len(i.split(";")[5].split(":")) > 1 else "UNK"
    )
)
akshay["akshay_motif_size"] = akshay["Info"].apply(
    lambda i: (
        len(i.split(";")[5].split(":")[1])
        if len(i.split(";")[5].split(":")) > 1
        else "UNK"
    )
)
akshay = akshay[akshay["akshay_motif"] != "UNK"]

# get mutations with a single motif, these are easy
single_motif_mutations = mutations[(mutations["n_motifs"] == 1) | (mutations["min_motiflen"] == mutations["max_motiflen"])]
single_motif_mutations["motif_sequence"] = single_motif_mutations["motifs"]
single_motif_mutations["consensus_motif_size"] = single_motif_mutations["min_motiflen"]

# merge multi-motif mutations with akshays file to decompose
multi_motif_mutations = mutations[mutations["n_motifs"] > 1]
multi_motif_mutations = akshay.merge(
    multi_motif_mutations,
    how="outer",
    left_on="ID",
    right_on="trid",
    indicator=True,
)

# get the multi motif mutations that were decomposable
decomposable_multi_motif = multi_motif_mutations[multi_motif_mutations["_merge"] == "both"].rename(
    columns={
        "akshay_motif": "motif_sequence",
        "akshay_motif_size": "consensus_motif_size",
    }
)
non_decomposable_multi_motif = multi_motif_mutations[
    multi_motif_mutations["_merge"] == "right_only"
]
non_decomposable_multi_motif["motif_sequence"] = non_decomposable_multi_motif["motifs"]
non_decomposable_multi_motif["consensus_motif_size"] = non_decomposable_multi_motif["min_motiflen"]

decomposable_multi_motif["decomposable"] = True
non_decomposable_multi_motif["decomposable"] = False

single_motif_mutations["decomposable"] = True

combined = pd.concat([single_motif_mutations, decomposable_multi_motif, non_decomposable_multi_motif])



mutations["consensus_motif_size"] = mutations["min_motiflen"]
mutations["motif_size_str"] = mutations["consensus_motif_size"].apply(lambda m: str(int(m)) if m < 7 else "7+")
mutations["decomposable"] = True



f, ax = plt.subplots(figsize=(3, 6))
counts = mutations.groupby("n_motifs").size().reset_index().rename(columns={0: "count"})
ax.bar(counts["n_motifs"], counts["count"], 1, lw=1, ec="w", color="cornflowerblue")
sns.despine(ax=ax)
ax.set_xlabel("Number of motifs\nin TR locus")
ax.set_ylabel("Number of loci")
ax.set_xticks(range(1, 6))
f.tight_layout()
f.savefig("motif_dist.png", dpi=200)


counts = mutations.groupby(["sample_id", "decomposable", "consensus_motif_size", "motif_size_str"]).size().reset_index().rename(columns={0: "count"})

counts["generation"] = counts["sample_id"].apply(
    lambda s: (
        "G4A"
        if str(s).startswith("2000")
        else "G4B" if str(s).startswith("2001") else "G3"
    ),
)

smps = counts["sample_id"].unique()
smp2idx = dict(zip(smps, range(len(smps))))

ind = np.arange(counts["sample_id"].nunique())

palette = ['006BA4', 'FF800E', 'ABABAB', '595959', '5F9ED1', 'C85200', '898989', 'A2C8EC', 'FFBC79', 'CFCFCF']
palette = [f"#{c}" for c in palette]
palette = sns.color_palette("colorblind", 7)

f, ax = plt.subplots(figsize=(7, 7))

bottom = np.zeros_like(ind)

# ax = axarr[dec_i]
for motif_i, (motif, motif_df) in enumerate(counts.groupby("motif_size_str")):
    motif_counts = np.zeros_like(ind)
    for s, i in smp2idx.items():
        sc = motif_df[motif_df["sample_id"] == s]["count"]
        if sc.shape[0] == 0: continue
        sc = sc.values[0]
        motif_counts[i] = sc
    ax.bar(
        ind,
        motif_counts,
        0.85,
        ec="w",
        lw=1,
        color=palette[motif_i],
        bottom=bottom,
        label=motif,
    )
    bottom += motif_counts
ax.set_xticks(ind)
ax.set_xticklabels(counts["sample_id"].unique(), rotation=45)
ax.set_xlabel("Sample ID")
ax.set_ylabel("# of DNMs")
ax.legend(title="Minimum motif size in locus (bp)", shadow=True, fontsize=10)

sns.despine(ax=ax)
f.tight_layout()
f.savefig("counts.png", dpi=200)
