import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as ss
from decimal import Decimal

plt.rc("font", size=12)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

def get_motif_types(row: pd.Series):
    if row["max_motiflen"] == 1:
        return "homopolymer"
    else:
        if row["simple_motif_size"] == "STR":
            return "non-homopolymer STR"
        else:
            return row["simple_motif_size"]
        
        
def is_parent_het(row: pd.Series, with_dnm: bool = True):
    poi = row["phase"]
    if poi == "dad":
        col = "father_AL" if with_dnm else "mother_AL"
        return False if len(set(row[col].split(","))) == 1 else True
    elif poi == "mom":
        col = "mother_AL" if with_dnm else "father_AL"
        return False if len(set(row[col].split(","))) == 1 else True


rng = np.random.default_rng(42)

ASSEMBLY = "CHM13v2"

mutations = pd.read_csv(f"{ASSEMBLY}.filtered.tsv", sep="\t", dtype={"sample_id": str})

mutations = mutations[(mutations["phase"] != "unknown") & (~mutations["#chrom"].isin(["chrX", "chrY"]))]


mutations["TR type"] = mutations.apply(lambda row: get_motif_types(row), axis=1)

# annotate every phased mutation with two columns. the first (`dnm`) will return
# False if the parent-of-origin is homozygous for TR AL, and True if the parent-of-origin
# is heterozygous. the second column (`wt`) will return False if the other parent (the
# parent who donated the non-denovo allele) is homozygous, and True if they're heterozygous.
mutations["transmitting_parent_is_het"] = mutations.apply(lambda row: is_parent_het(row, with_dnm=True), axis=1)
mutations["non_transmitting_parent_is_het"] = mutations.apply(lambda row: is_parent_het(row, with_dnm=False), axis=1)

tidy = mutations[
    [
        "trid",
        "sample_id",
        "TR type",
        "transmitting_parent_is_het",
        "non_transmitting_parent_is_het",
    ]
].melt(
    id_vars=["trid", "TR type"],
    value_vars=["transmitting_parent_is_het", "non_transmitting_parent_is_het"],
)

tidy = tidy.groupby(["TR type", "variable"]).agg(count = ("value", "sum")).reset_index()
motif_counts = mutations.groupby("TR type").size().reset_index().rename(columns={0: "total"})

tidy = tidy.merge(motif_counts)

tidy["Fraction of DNMs"] = tidy["count"] / tidy["total"]

tidy["non_count"] = tidy["total"] - tidy["count"] 

# tidy["variable"] = tidy["variable"].apply(lambda v: "Transmitted DNM" if v == "dnm" else "Did not transmit DNM")

f, axarr = plt.subplots(1, 4, figsize=(9, 7), sharex=True, sharey=False)

tr2idx = dict(zip(tidy["TR type"].unique(), range(tidy["TR type"].nunique())))

colors = sns.color_palette("colorblind")

for mi, (motif, motif_df) in enumerate(tidy.groupby("TR type")):
    motif_df = motif_df.sort_values("variable", ascending=False)
    ax = axarr[tr2idx[motif]]

    contingency = motif_df[["count", "non_count"]].values
    
    res = ss.chi2_contingency(contingency)
    pval = "{:.1e}".format(Decimal(str(res.pvalue))) if res.pvalue != 1 else "1"

    vals = motif_df["Fraction of DNMs"]
    top = 1 - vals
    ax.bar(motif_df["variable"], vals, ec="w", lw=2, label="Heterozygous", color=colors[0])
    ax.bar(motif_df["variable"], top, bottom=vals, ec="w", lw=2, label="Homozygous", color="gainsboro")

    for vi, v in enumerate(vals):
        ax.text(vi - 0.1, v / 2, f"{contingency[vi, 0]}", family="monospace", color="w")
        ax.text(vi - 0.1, v + ((1 - v) / 2.2), f"{contingency[vi, 1]}", family="monospace", color="k")
    ax.set_xticklabels(["DN", "NDN"])
    ax.set_title(f"{motif}\n" + r"$\chi^{2}$" + f" p = {pval}")
    sns.despine(ax=ax, left=mi > 0)
    if mi > 0:
        ax.set_yticks([])
    if mi == 0:
        ax.set_ylabel("Fraction of loci")
    if mi == 1:
        ax.set_xlabel("Which allele did this parent transmit?")
    
axarr[0].legend(title="Parental genotype", shadow=True)
f.tight_layout()
f.savefig("heterozygote_effect.png", dpi=200)
