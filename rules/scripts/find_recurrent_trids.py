import pandas as pd
from collections import Counter
from cyvcf2 import VCF
import numpy as np


from snakemake.script import snakemake

ORIG = """chr1_54393726_54394070_trsolve
chr8_2376919_2377075_trsolve
chr7_2500010_2500042_trsolve
chr4_79949242_79949442_trsolve
chr4_21696993_21697153_trsolve
chr12_119907035_119907158_trsolve
chr12_114852499_114852706_trsolve
chr7_42892201_42892385_trsolve
chr21_33731357_33731465_trsolve
chr9_36529968_36530006_trsolve
chr7_6540708_6540973_trsolve
chr7_152489617_152489683_trsolve
chr12_95884953_95885246_trsolve
chr7_13334154_13334671_trsolve
chr15_32243116_32243499_trsolve
chr14_95031468_95031513_trsolve
chrX_49532133_49532368_trsolve
chrX_45925532_45926294_trsolve
chr8_91356817_91357059_trsolve
chr7_2490933_2492021_trsolve
chr7_19629376_19630006_trsolve
chr6_51818304_51818650_trsolve
chr2_121116093_121116213_trsolve
chr2_102413803_102413989_trsolve
chr20_38333021_38333121_trsolve
chr1_153340368_153340485_trsolve
chr19_6620632_6620703_trsolve
chr19_23902332_23902399_trsolve
chr13_71023191_71023231_trsolve
chr11_56923470_56923565_trsolve
chr10_409491_410951_trsolve
chr10_2390774_2391407_trsolve""".split()


dfs = []
for fh in snakemake.input.fhs:
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    dfs.append(df)
dfs = pd.concat(dfs)


samples_with_denovo = dfs.groupby("trid").agg(samples_with_denovo=("sample_id", lambda s: ",".join(s))).reset_index()
dfs = dfs.merge(samples_with_denovo)

# remove recurrents where one sample has multiple dnms
dfs["uniq_smp"] = dfs["samples_with_denovo"].apply(lambda s: len(set(s.split(","))))
# dfs = dfs[dfs["uniq_smp"] > 1]

dfs["generation"] = dfs["sample_id"].apply(lambda s: "G4" if s.startswith("200") else "G2" if s in ("2209", "2188") else "G3")

# figure out which TRIDs are observed as DNs
# multiple times in a single generation, rather than multiple
# times across generations
generational = (
    dfs.groupby("trid")
    .agg(
        inter=("generation", lambda g: len(set(g))),
        intra=("generation", lambda g: Counter(g).most_common()[0][1]),
    )
    .reset_index()
)

intra_trids = generational[generational["intra"] > 1]["trid"].unique()
inter_trids = generational[generational["inter"] > 1]["trid"].unique()
combined_trids = set(intra_trids).union(set(inter_trids))

print ("INTRA-GENERATIONAL", len(intra_trids))
print ("INTER-GENERATIONAL", len(inter_trids))
print ("UNIQUE", len(combined_trids))

# filter to DNMs that are recurrent across generations
dfs = dfs[dfs["trid"].isin(combined_trids)].drop_duplicates("trid")
dfs = dfs.merge(generational)


# for each recurrent DNM, make sure we have sufficient read depth in every member of the pedigree
vcfs = [VCF(fh, gts012=True) for fh in snakemake.input.vcfs]
res = []
for i, row in dfs.iterrows():
    chrom, start, end = row["#chrom"], row["start"], row["end"]
    region = f"{chrom}:{start}-{end}"
    trid = row["trid"]
    n_suff_dp = 0
    for vcf in vcfs:
        for v in vcf(region):
            if v.INFO.get("TRID") != trid: continue
            assert v.INFO.get("TRID") == trid
            gts = v.gt_types
            dp = np.sum(v.format("SD"))
            assert gts.shape[0] == 1
            
            if gts[0] == 3: continue
            if dp < 10: continue
            n_suff_dp += 1
    row_dict = row.to_dict()
    if n_suff_dp < len(vcfs):
        row_dict.update({"sufficient_cohort_depth": 0})
    else:
        row_dict.update({"sufficient_cohort_depth": 1})

    res.append(row_dict)

res_df = pd.DataFrame(res)

res_df[
    [
        "trid",
        "#chrom",
        "start",
        "end",
        "index",
        "motifs",
        "child_AL",
        "samples_with_denovo",
        "sufficient_cohort_depth",
        "inter",
        "intra",
    ]
].drop_duplicates("trid").to_csv(snakemake.output.out, index=False, sep="\t")
