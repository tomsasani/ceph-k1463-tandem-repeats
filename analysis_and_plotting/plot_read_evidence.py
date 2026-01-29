import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
from bx.intervals.intersection import Interval, IntervalTree
from collections import defaultdict
import csv

# get mutations
mutations = []
for fh in snakemake.input.mutations:
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    mutations.append(df)
mutations = pd.concat(mutations).dropna(subset=["kid_evidence", "mom_evidence", "dad_evidence"])


recurrents = pd.read_csv(snakemake.input.recurrent, sep="\t")["trid"].to_list()
mutations = mutations[~mutations["trid"].isin(recurrents)]

mutations["reference_al"] = mutations["end"] - mutations["start"]

mutations["denovo_al"] = mutations.apply(lambda row: list(map(int, row["child_AL"].split(",")))[row["index"]], axis=1)

# calculate child-to-reference diff
mutations["child_to_reference_diff"] = (
    mutations["denovo_al"] - mutations["reference_al"]
)
mutations["child_to_parent_diff"] = mutations["likely_denovo_size"]

mutations = mutations[mutations["min_motiflen"] >= snakemake.params.minimum_motif_size]
mutations = mutations[mutations["likely_denovo_size"] >= snakemake.params.minimum_dnm_size]


for (sample, trid, genotype), trid_df in mutations.groupby(
    ["alt_sample_id", "trid", "genotype"]
):
    validation = trid_df["validation_status"].unique()[0]
    motif_size = trid_df["min_motiflen"].values[0]
    poi, _, poi_support = trid_df["phase_consensus"].unique()[0].split(":")
    size = trid_df["likely_denovo_size"].unique()[0]

    # if motif_size > 100:
    #     print (trid_df)

    if poi != "unknown":
        poi = poi if float(poi_support) > 0.75 else "unknown"

    hoi, _, hoi_support = trid_df["haplotype_in_parent_consensus"].unique()[0].split(":")
    if hoi != "unknown":
        hoi = hoi if float(hoi_support) > 0.75 else "unknown"

    precursor_allele_length = trid_df["precursor_allele_length_in_parent"].unique()[0]
    ref_len = trid_df["end"].unique()[0] - trid_df["start"].unique()[0]
    precursor_allele_length = precursor_allele_length - ref_len
    chrom, start, end, _ = trid.split("_")

    is_male_sex_chrom = "," not in trid_df["child_AL"].unique()[0]
    
    if not is_male_sex_chrom:
        trid_df["non_denovo_al"] = trid_df.apply(
            lambda row: row["child_AL"].split(",")[1 - row["index"]], axis=1
        )
    if trid_df.shape[0] != 1: 
        continue

    mom_diffs, dad_diffs, kid_diffs = [], [], []
    for col in ("mom_evidence", "dad_evidence", "kid_evidence"):
        for diff_count in trid_df[col].values[0].split("|"):
            diff, count = list(map(int, diff_count.split(":")))
            if col == "mom_evidence":
                mom_diffs.extend([diff] * count)
            elif col == "dad_evidence":
                dad_diffs.extend([diff] * count)
            else:
                kid_diffs.extend([diff] * count)

    exp_denovo_diff = trid_df["exp_allele_diff_denovo"].unique()[0]
    exp_non_denovo_diff = trid_df["exp_allele_diff_non_denovo"].unique()[0]

    exp_denovo = trid_df["denovo_al"].unique()[0]
    if not is_male_sex_chrom:
        exp_non_denovo = trid_df["non_denovo_al"].unique()[0]

    f, ax = plt.subplots()

    plot_df = []
    for diff_list, label in zip(
        (kid_diffs, mom_diffs, dad_diffs), ("kid", "mom", "dad")
    ):
        is_poi = None
        if label == "kid":
            is_poi = "kid"
        else:
            if poi == "unknown":
                is_poi = False
            else:
                is_poi = poi == label
        for diff in diff_list:
            plot_df.append({"sample": label, "is_poi": is_poi, "diff": diff})
    plot_df = pd.DataFrame(plot_df)
    sns.stripplot(
        data=plot_df,
        x="sample",
        y="diff",
        ax=ax,
        hue="is_poi",
        palette={True: "green", False: "red", "kid": "blue"},
    ) 
    ax.axhline(
        exp_denovo_diff,
        ls=":",
        c="dodgerblue",
        label="Expected allele length (de novo)",
    )
    if not is_male_sex_chrom:
        ax.axhline(
            exp_non_denovo_diff,
            ls=":",
            c="firebrick",
            label="Expected allele length (non de novo)",
        )
    if precursor_allele_length != -1:
        ax.axhline(precursor_allele_length, ls="--", c="seagreen", label="Expected precursor allele in parent")
    ax.set_ylabel(f"Allele length inferred using {snakemake.wildcards.TECH} reads")
    ax.set_xlabel("Sample")
    ax.set_xticks(range(3))
    ax.set_xticklabels(["Kid", "Mom", "Dad"])
    ax.set_title(f"{trid}\nMotif size = {motif_size}\nDNM size = {size}")
    if abs(exp_non_denovo_diff) > 500 or abs(exp_denovo_diff) > 500:
        ax.set_yscale("symlog")
    # ax.legend()
    sns.despine(ax=ax)
    f.tight_layout()
    f.savefig(
        f"{snakemake.params.outpref}/{snakemake.wildcards.TECH}.{snakemake.wildcards.ASSEMBLY}.{motif_size}.{sample}.{trid}.png",
        dpi=200,
    )
    plt.close()

with open(snakemake.output.fh, "w") as outfh:
    for l in trid_df["trid"].unique():
        print (l, file=outfh)
outfh.close()
