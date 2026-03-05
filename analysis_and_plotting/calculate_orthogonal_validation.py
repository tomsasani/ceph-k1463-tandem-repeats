import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from bx.intervals.intersection import Interval, IntervalTree
import csv
from snakemake.script import snakemake

def annotate_with_censat(row: pd.Series, censat):
    overlaps = censat[row["#chrom"]].find(row["start"], row["end"])
    if len(overlaps) == 0: 
        return "no"
    else:
        return overlaps[0].value["kind"].split("_")[0]


censat = defaultdict(IntervalTree)
with open(snakemake.input.censat, "r") as infh:
    csvf = csv.reader(infh, delimiter="\t")
    for l in csvf:
        chrom, start, end = l[:3]
        censat[chrom].insert_interval(Interval(int(start), int(end), value={"kind": l[3]}))


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

mutations = mutations[mutations["max_motiflen"] > 1]

if snakemake.wildcards.ASSEMBLY == "CHM13v2":
    mutations["overlaps_censat"] = mutations.apply(lambda row: annotate_with_censat(row, censat), axis=1)
else:
    mutations["overlaps_censat"] = "no"


mutations["is_phased"] = mutations["phase_consensus"].apply(lambda p: "Y" if "unknown" not in p else "N")
mutations["is_phased_int"] = mutations["is_phased"].apply(lambda p: 1 if p == "Y" else 0)

mutations = mutations[mutations["validation_status"] != "no_data"]

mutations = mutations[~mutations["trid"].isin(recurrent_trids)]

mutations["in_censat"] = mutations["overlaps_censat"].apply(lambda s: "N" if s == "no" else "Y")

mutations["is_validated"] = mutations["validation_status"].apply(lambda v: 1 if v == "pass" else 0)

group_cols = ["simple_motif_size", "in_censat", "validation_status"]


counts = mutations.groupby(group_cols).agg(passing=("is_validated", "sum"), total=("is_validated", "count")).reset_index()
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

f, ax = plt.subplots(figsize=(5, 8))
sns.barplot(
    data=res,
    y="frac",
    x="simple_motif_size",
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
