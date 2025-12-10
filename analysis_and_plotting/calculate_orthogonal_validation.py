import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


plt.rc("font", size=14)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

recurrent = pd.read_csv(snakemake.input.recurrent, sep="\t")
recurrent_trids = recurrent["trid"].unique()

# get mutations
mutations = []
for fh in snakemake.input.mutations:
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    mutations.append(df)
mutations = pd.concat(mutations)
# mutations = mutations[mutations["simple_motif_size"].isin(["homopolymer", "non-homopolymer STR"]]

print (mutations.groupby(["simple_motif_size", "validation_status"]).size())
mutations["is_phased"] = mutations["phase_consensus"].apply(lambda p: "Y" if "unknown" not in p else "N")
mutations["is_phased_int"] = mutations["is_phased"].apply(lambda p: 1 if p == "Y" else 0)

mutations = mutations[mutations["validation_status"] != "no_data"]

print (mutations.shape)
mutations = mutations[~mutations["trid"].isin(recurrent_trids)]
print (mutations.shape)

mutations["TR type"] = mutations["simple_motif_size"]
mutations["is_expansion"] = mutations["likely_denovo_size"].apply(lambda s: "expansion" if s >= 1 else "contraction")

mutations["is_validated"] = mutations["validation_status"].apply(lambda v: 1 if v == "pass" else 0)

group_cols = ["TR type", "is_expansion"]

counts = mutations.groupby(group_cols).agg(passing=("is_validated", "sum")).reset_index()
chi2 = [list(counts.values[0:2, 2]), list(counts.values[2:, 2])]


bootstraps = 1_000
res = []
for bs in range(bootstraps):
    _mutations = mutations.sample(frac=1, replace=True)
    counts = _mutations.groupby(group_cols).agg(passing=("is_validated", "sum")).reset_index()
    totals = _mutations.groupby(group_cols).size().reset_index().rename(columns={0: "total"})
    counts = counts.merge(totals)
    counts["trial"] = bs
    res.append(counts)

res = pd.concat(res)
res["frac"] = res["passing"] / res["total"]

res = res.rename(columns={"is_expansion": "Size"})

f, ax = plt.subplots(figsize=(5, 8))
sns.barplot(
    data=res,
    hue="Size",
    y="frac",
    x="TR type",
    estimator="mean",
    errorbar=("ci", 95),
    ec="w",
    lw=1,
    ax=ax,
)
ax.set_ylabel("Element validation rate")
sns.despine(ax=ax)
f.tight_layout()
f.savefig(snakemake.output.png, dpi=200)
