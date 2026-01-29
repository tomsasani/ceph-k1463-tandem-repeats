import pandas as pd
from bx.intervals.intersection import Interval, IntervalTree
from collections import defaultdict
import csv
from analysis_and_plotting.FAILING_TRIDS import FAIL_VNTRS, FAIL_SVS


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

# read in all per-sample DNM files
mutations = []
for fh in snakemake.input.dnms:
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    mutations.append(df)

mutations = pd.concat(mutations).fillna(
    {
        "children_with_denovo_allele": "unknown",
        "grandparents_with_denovo_allele": "unknown",
    }
)

mutations = mutations[~mutations["trid"].isin(FAIL_VNTRS)]
mutations = mutations[~mutations["trid"].isin(FAIL_SVS)]

if snakemake.wildcards.ASSEMBLY == "CHM13v2":
    mutations["overlaps_censat"] = mutations.apply(lambda row: annotate_with_censat(row, censat), axis=1)
else:
    mutations["overlaps_censat"] = "no"

# get sample IDs so we can filter the denominator files
sample_ids = mutations["sample_id"].unique()
# map alternate (NAXXXX) IDs to original (2189) IDs
alt2orig = dict(zip(mutations["alt_sample_id"], mutations["sample_id"]))
orig2alt = {v:k for k,v in alt2orig.items()}

# if desired, remove recurrent sites that are recurrent in G3
if snakemake.params.filter_by_recurrent:
    recurrent = pd.read_csv(snakemake.input.recurrent, sep="\t")
    recurrent_trids = recurrent["trid"].unique()
    mutations = mutations[~mutations["trid"].isin(recurrent_trids)]

# if desired, filter DNMs that didn't pass our Element validation checks.
# element validation data should already be in this dataframe as an extra annotation.
if snakemake.params.filter_by_orthogonal:
    # filter to STRs that had orthogonal data
    orthogonal = mutations[
        (mutations["simple_motif_size"].isin(["STR", "non-homopolymer STR", "homopolymer"]))
        & (mutations["validation_status"] != "no_data")
    ]
    # make sure we're only looking at sites with max AL <= 120
    orthogonal["max_al"] = orthogonal["child_AL"].apply(lambda a: max(map(int, a.split(","))))
    orthogonal = orthogonal[orthogonal["max_al"] <= 120]

    fail_trids = orthogonal[orthogonal["validation_status"] != "pass"]["trid"].unique()

    mutations = mutations[~mutations["trid"].isin(fail_trids)]

# if desired, remove untransmitted DNMs in the samples for whom we can assess that
if snakemake.params.filter_by_transmission:

    has_transmission = mutations[mutations["sample_id"].isin(["2189", "2216"])]
    is_transmitted = has_transmission[has_transmission["children_with_denovo_allele"] != "unknown"]
    mutations = mutations[
        (~mutations["sample_id"].isin(["2189", "2216"]))
        | (mutations["children_with_denovo_allele"] != "unknown")
    ]

if snakemake.params.filter_by_grandparents:
    mutations = mutations[
        (mutations["sample_id"].str.startswith("200"))
        | (mutations["grandparents_with_denovo_allele"] == "unknown")
    ]

mutations["generation"] = mutations["paternal_id"].apply(
    lambda s: (
        "G4A"
        if s == 200080
        else (
            "G4B"
            if s == 2189
            else "G3" if s == 2209 else "G2A" if s == 2281 else "G2B"
        )
    )
)

# output filtered mutations to CSV
mutations.to_csv(snakemake.output.csv, index=False, sep="\t")
