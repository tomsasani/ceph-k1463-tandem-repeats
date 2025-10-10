import pandas as pd

include: "rules/alignment.smk"
include: "rules/snvs.smk"
include: "rules/trgt.smk"
include: "rules/phasing.smk"
include: "rules/trgt_denovo.smk"
include: "rules/denovo_mutation_filtering.smk"


# create global dictionaries we'll use
PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str})


CHILDREN = ped[(ped["paternal_id"] != "missing") & (~ped["paternal_id"].isin(["2281", "2214"]))]["sample_id"].to_list()

rule all:
    input:
        expand("csv/filtered_and_merged/{SAMPLE}.CHM13v2.TOPUP.tsv", SAMPLE=["200101"]),
        expand("csv/annotated/{SAMPLE}.CHM13v2.element.TOPUP.tsv", SAMPLE=["200101"]),
        expand("csv/recurrent/CHM13v2.TOPUP.recurrent.tsv")
