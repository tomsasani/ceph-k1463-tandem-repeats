import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from bx.intervals.intersection import Interval, IntervalTree
from collections import defaultdict
import csv


def annotate_with_censat(row: pd.Series, censat):
    overlaps = censat[row["#chrom"]].find(row["start"], row["end"])
    if len(overlaps) == 0: 
        return "no"
    else:
        return overlaps[0].value["kind"].split("_")[0]


censat = defaultdict(IntervalTree)
with open("data/t2t.censat.bed", "r") as infh:
    csvf = csv.reader(infh, delimiter="\t")
    for l in csvf:
        chrom, start, end = l[:3]
        censat[chrom].insert_interval(Interval(int(start), int(end), value={"kind": l[3]}))


pd.set_option("display.precision", 8)
plt.rc("font", size=16)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

# define some global variables we'll access for filtering
FILTER_TRANSMITTED = False
FILTER_GRANDPARENTS = False
PER_HAPLOTYPE = True
FILTER_RECURRENT = True
FILTER_ELEMENT = True

# define the assembly we're sourcing our DNMs from
ASSEMBLY = "GRCh38"
TOPUP = "TOPUP"
ORTHOGONAL_TECH = "element"

# read in all per-sample DNM files
mutations = []
for fh in glob.glob(f"csv/annotated/*.{ASSEMBLY}.{ORTHOGONAL_TECH}.{TOPUP}.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    mutations.append(df)

mutations = pd.concat(mutations).fillna(
    {
        "children_with_denovo_allele": "unknown",
        "grandparents_with_denovo_allele": "unknown",
    }
)

if ASSEMBLY == "CHM13v2":
    mutations["overlaps_censat"] = mutations.apply(lambda row: annotate_with_censat(row, censat), axis=1)
else:
    mutations["overlaps_censat"] = 0


# get sample IDs so we can filter the denominator files
sample_ids = mutations["sample_id"].unique()
# map alternate (NAXXXX) IDs to original (2189) IDs
alt2orig = dict(zip(mutations["alt_sample_id"], mutations["sample_id"]))
orig2alt = {v:k for k,v in alt2orig.items()}

# read in per-sample denominators
denoms = []
for fh in glob.glob(f"csv/denominators/*.{ASSEMBLY}.v4.0.denominator.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})    
    denoms.append(df)
denoms = pd.concat(denoms)
denoms = denoms[denoms["sample_id"].isin(sample_ids)]

# if we want to calculate mutation rates per-haplotype, rather
# than per-genome, we need to multiply the denominator by 2 to account
# for there being 2 copies of every locus in a diploid genome.
if PER_HAPLOTYPE:
    denoms["denominator"] *= 2

# if desired, remove recurrent sites that are recurrent in G3
if FILTER_RECURRENT:
    recurrent = pd.read_csv(f"csv/recurrent/{ASSEMBLY}.v4.0.recurrent.tsv", sep="\t")
    # recurrent["contains_G3"] = recurrent["samples_with_denovo"].apply(lambda samples: any([s in sample_ids for s in samples.split(",")]))
    # recurrent = recurrent[recurrent["contains_G3"] == True]
    recurrent_trids = recurrent["trid"].unique()
    mutations = mutations[~mutations["trid"].isin(recurrent_trids)]

# if desired, filter DNMs that didn't pass our Element validation checks.
# element validation data should already be in this dataframe as an extra annotation.
if FILTER_ELEMENT:
    # filter to STRs that had orthogonal data
    orthogonal = mutations[
        (mutations["simple_motif_size"] == "STR")
        & (mutations["validation_status"] != "no_data")
    ]
    # make sure we're only looking at sites with max AL <= 120
    orthogonal["max_al"] = orthogonal["child_AL"].apply(lambda a: max(map(int, a.split(","))))
    orthogonal = orthogonal[orthogonal["max_al"] <= 120]

    fail_trids = orthogonal[orthogonal["validation_status"] != "pass"]["trid"].unique()

    mutations = mutations[~mutations["trid"].isin(fail_trids)]

# if desired, remove untransmitted DNMs in the samples for whom we can assess that
if FILTER_TRANSMITTED:

    has_transmission = mutations[mutations["sample_id"].isin(["2189", "2216"])]

    is_transmitted = has_transmission[has_transmission["children_with_denovo_allele"] != "unknown"]

    mutations = mutations[
        (~mutations["sample_id"].isin(["2189", "2216"]))
        | (mutations["children_with_denovo_allele"] != "unknown")
    ]

if FILTER_GRANDPARENTS:
    mutations = mutations[
        (mutations["sample_id"].str.startswith("200"))
        | (mutations["grandparents_with_denovo_allele"] == "unknown")
    ]

phased = mutations[mutations["phase"] != "unknown"]

hap_phased = phased[~phased["haplotype_in_parent_consensus"].str.contains("unknown")]

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


# figure out which DNMs occurred at complex loci
complex_mutations = mutations[mutations["n_motifs"] > 1]
simple_mutations = mutations[mutations["n_motifs"] == 1]

# read in akshay file
akshay = pd.read_csv("data/K1463.CHM13v2.DNMs.416.demintr.output", sep="\t")

# merge complex loci with akshay's file
complex_decomposed = complex_mutations.merge(akshay, left_on="trid", right_on="ID")

# figure out which of these sites have only a single mutational event
n_events = complex_decomposed.groupby("ID").size().reset_index().rename(columns={0: "n_events"})
complex_decomposed = complex_decomposed.merge(n_events)

def get_mutating_motif(info):
    parent_repeat = info.split(";")[5]
    parent_repeat_seq = parent_repeat.split(":")
    if len(parent_repeat_seq) > 1:
        return parent_repeat.split(":")[1]
    else:
        return "UNK"

complex_decomposed = complex_decomposed[complex_decomposed["n_events"] == 1]

complex_decomposed["akshay_motif"] = complex_decomposed["Info"].apply(lambda i: get_mutating_motif(i))
complex_decomposed["akshay_motif_size"] = complex_decomposed["akshay_motif"].apply(lambda m: len(m) if m != "UNK" else m)
# mutations = mutations[mutations["precursor_sequence_in_parent"] == "unknown"]
# mutations["parent_of_origin"] = mutations["phase_consensus"].apply(lambda p: p.split(":")[0])


# figure out which loci only contain one kind of motif -- regardless of whether we
# are able to determine a haplotype of origin, we can still know the exact size
# of the motif that mutated.


# output filtered mutations to CSV
mutations.to_csv(f"{ASSEMBLY}.filtered.tsv", index=False, sep="\t")
