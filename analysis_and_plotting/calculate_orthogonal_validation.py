import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import numpy as np
import scipy.stats as ss

from bx.intervals.intersection import Interval, IntervalTree
from collections import defaultdict
import csv

def annotate_with_censat(row: pd.Series, censat):
    overlaps = censat[row["#chrom"]].find(row["start"], row["end"])
    if len(overlaps) == 0: 
        return "no"
    else:
        return overlaps[0].value["kind"].split("_")[0]

plt.rc("font", size=14)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]


TECH = "element"
ASSEMBLY = "GRCh38"
MIN_SIZE = 1
FILTER_RECURRENT = False
USE_NEW_BAM = "TOPUP"

recurrent = pd.read_csv(f"csv/recurrent/{ASSEMBLY}.{USE_NEW_BAM}.recurrent.tsv", sep="\t")
recurrent_trids = recurrent["trid"].unique()

censat = defaultdict(IntervalTree)
with open("data/t2t.censat.bed", "r") as infh:
    csvf = csv.reader(infh, delimiter="\t")
    for l in csvf:
        chrom, start, end = l[:3]
        censat[chrom].insert_interval(Interval(int(start), int(end), value={"kind": l[3]}))

# get mutations
mutations = []
for fh in glob.glob(f"csv/annotated/*.{ASSEMBLY}.{TECH}.{USE_NEW_BAM}.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    mutations.append(df)
mutations = pd.concat(mutations)
mutations = mutations[mutations["simple_motif_size"] == "STR"]

if ASSEMBLY == "CHM13v2":
    mutations["overlaps_censat"] = mutations.apply(lambda row: annotate_with_censat(row, censat), axis=1)
else:
    mutations["overlaps_censat"] = 0

if FILTER_RECURRENT:
    mutations = mutations[~mutations["trid"].isin(recurrent_trids)]


mutations["is_phased"] = mutations["phase_consensus"].apply(lambda p: "Y" if "unknown" not in p else "N")
mutations["is_phased_int"] = mutations["is_phased"].apply(lambda p: 1 if p == "Y" else 0)

mutations = mutations[mutations["validation_status"] != "no_data"]

mutations["TR type"] = mutations.apply(
    lambda row: (
        "Homopolymer"
        if 1 in (row["min_motiflen"], row["max_motiflen"])
        else "non-homopolymer STR"
    ),
    axis=1,
)
mutations["is_expansion"] = mutations["likely_denovo_size"].apply(lambda s: "expansion" if s >= 1 else "contraction")

# mutations["kid_tot"] = mutations["kid_evidence"].apply(lambda e: sum([int(_e.split(":")[1]) for _e in e.split("|")]))
# mutations["mom_tot"] = mutations["mom_evidence"].apply(lambda e: sum([int(_e.split(":")[1]) for _e in e.split("|")]))
# mutations["dad_tot"] = mutations["dad_evidence"].apply(lambda e: sum([int(_e.split(":")[1]) for _e in e.split("|")]))

mutations["is_validated"] = mutations["validation_status"].apply(lambda v: 1 if v == "pass" else 0)

print (mutations.shape)
group_cols = ["TR type", "is_expansion"]

counts = mutations.groupby(group_cols).agg(passing=("is_validated", "sum")).reset_index()
chi2 = [list(counts.values[0:2, 2]), list(counts.values[2:, 2])]
print (ss.chi2_contingency(chi2))
bootstraps = 100
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
f.savefig("ortho.png", dpi=200)
