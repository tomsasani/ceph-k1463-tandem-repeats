import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as ss
from decimal import Decimal

plt.rc("font", size=13)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]


        
def is_parent_het(row: pd.Series, with_dnm: bool = True):
    poi = row["phase"]
    if poi == "dad":
        col = "father_AL" if with_dnm else "mother_AL"
        return False if len(set(row[col].split(","))) == 1 else True
    elif poi == "mom":
        col = "mother_AL" if with_dnm else "father_AL"
        return False if len(set(row[col].split(","))) == 1 else True

def allele_length_diff(row: pd.Series, with_dnm: bool = True):
    poi = row["phase"]
    diff = -1
    if poi == "dad":
        col = "father_AL" if with_dnm else "mother_AL"
        als = list(map(int, row[col].split(",")))
        diff = abs(als[1] - als[0])
    elif poi == "mom":
        col = "mother_AL" if with_dnm else "father_AL"
        als = list(map(int, row[col].split(",")))
        diff = abs(als[1] - als[0])
    return diff


rng = np.random.default_rng(42)

mutations = pd.read_csv(snakemake.input.mutations, sep="\t", dtype={"sample_id": str})

mutations = mutations[(mutations["phase"] != "unknown") & (~mutations["#chrom"].isin(["chrX", "chrY"]))]


# NOTE: for these plots, we discriminate between homopolymer and non-homopolymer
# containing loci, because many prior studies have focused on the former.
def adjusted_tr_type(row):
    if row["max_motiflen"] == 1:
        return "homopolymer"
    else:
        # this means that a locus with a homopolymer and another STR will be "STR"
        return row["simple_motif_size"]

mutations["TR type"] = mutations.apply(lambda row: adjusted_tr_type(row), axis=1)


# annotate every phased mutation with two columns. the first (`dnm`) will return
# False if the parent-of-origin is homozygous for TR AL, and True if the parent-of-origin
# is heterozygous. the second column (`wt`) will return False if the other parent (the
# parent who donated the non-denovo allele) is homozygous, and True if they're heterozygous.
mutations["transmitting_parent_is_het"] = mutations.apply(lambda row: is_parent_het(row, with_dnm=True), axis=1)
mutations["non_transmitting_parent_is_het"] = mutations.apply(lambda row: is_parent_het(row, with_dnm=False), axis=1)


mutations["diff_poi"] = mutations.apply(lambda row: allele_length_diff(row, with_dnm=True), axis=1)
mutations["diff_npoi"] = mutations.apply(lambda row: allele_length_diff(row, with_dnm=False), axis=1)

f, ax = plt.subplots()

mutations["diff"] = mutations["diff_npoi"] - mutations["diff_poi"]
max_diff = np.max(np.abs(mutations["diff"].values))
# mutations = mutations[mutations["TR type"] == "homopolymer"]
for i,row in mutations.iterrows():
    x1, x2 = 0, 1
    y1, y2 = row["diff_npoi"], row["diff_poi"]
    # jitter x
    x1 += np.random.normal(0, 0.01)
    x2 += np.random.normal(0, 0.01)
    ax.scatter([x1, x2], [y1, y2], c="dodgerblue")
    alpha = abs(y2 - y1) / max_diff
    ax.plot([x1, x2], [y1, y2], c="grey", alpha=alpha)
ax.set_xlim(-0.5, 1.5)
f.savefig("o.png")

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

f, axarr = plt.subplots(1, 4, figsize=(10, 7), sharex=True, sharey=False)

motifs = ["homopolymer", "STR", "VNTR", "STR + VNTR"]

tr2idx = dict(zip(motifs, range(len(motifs))))

colors = sns.color_palette("colorblind")

for motif, motif_df in tidy.groupby("TR type"):
    mi = tr2idx[motif]
    motif_df = motif_df.sort_values("variable", ascending=False)

    ax = axarr[tr2idx[motif]]

    contingency = motif_df[["count", "non_count"]].values
    
    res = ss.chi2_contingency(contingency)
    pval = "{:.1e}".format(Decimal(str(res.pvalue))) if res.pvalue != 1 else "1"

    vals = motif_df["Fraction of DNMs"]
    top = 1 - vals
    ax.bar(motif_df["variable"], vals, ec="w", color=["lightsteelblue", "gainsboro"], lw=2, label="Heterozygous")#, color="dimgray")
    ax.bar(motif_df["variable"], top, bottom=vals, ec="w", lw=2, label="Homozygous", color=["#3274a2", "dimgray"])

    for vi, v in enumerate(vals):
        v1, v2 = contingency[vi, 0], contingency[vi, 1]
        ax.text(vi - (0.05 * len(str(v1))), v / 2, f"{v1}", family="monospace", color="k")
        ax.text(vi - (0.05 * len(str(v2))), v + ((1 - v) / 2.2), f"{v2}", family="monospace", color="w")
    ax.set_xticklabels(["PO", "NPO"])
    ax.set_title(f"{motif}\n" + r"$\chi^{2}$" + f" p = {pval}")
    sns.despine(ax=ax, left=mi > 0)
    if mi > 0:
        ax.set_yticks([])
    if mi == 0:
        ax.set_ylabel("Fraction of loci")
    if mi == 1:
        ax.set_xlabel("")
    
axarr[0].legend(title="Parental genotype", shadow=True, fontsize=14, loc="lower left")
f.tight_layout()
f.savefig(snakemake.output.png, dpi=200)
