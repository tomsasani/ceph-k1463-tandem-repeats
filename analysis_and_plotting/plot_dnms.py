import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from snakemake.script import snakemake


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]
plt.rc("font", size=14)


mutations = pd.read_csv(snakemake.input.mutations, sep="\t", dtype={"sample_id": str, "paternal_id": str})

mutations["motif_size_str"] = mutations["min_motiflen"].apply(lambda m: str(int(m)) if m < 7 else "7+")
mutations["decomposable"] = True

counts = (
    mutations.groupby(
        [
            "sample_id",
            "motif_size_str",
        ]
    )
    .size()
    .reset_index()
    .rename(columns={0: "count"})
)

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

f, ax = plt.subplots(figsize=(9, 7))

bottom = np.zeros_like(ind)

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
ax.legend(title="Minimum motif size in locus (bp)", shadow=True, fontsize=12)

sns.despine(ax=ax)
f.tight_layout()
f.savefig(snakemake.output.png, dpi=200)
