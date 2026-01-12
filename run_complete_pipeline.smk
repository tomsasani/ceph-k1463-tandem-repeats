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
ALL_SAMPLES = ped["sample_id"].to_list()
ASSEMBLIES = ["CHM13v2"]


rule all:
    input:
        expand("csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.TOPUP.tsv", SAMPLE=CHILDREN, ASSEMBLY=ASSEMBLIES),
        expand("csv/denominators/{SAMPLE}.{ASSEMBLY}.TOPUP.denominator.tsv", SAMPLE=CHILDREN, ASSEMBLY=ASSEMBLIES),
        expand("csv/annotated/{SAMPLE}.{ASSEMBLY}.{TECH}.TOPUP.tsv", SAMPLE=CHILDREN, ASSEMBLY=ASSEMBLIES, TECH=["element", "hifi"]),
        expand("csv/recurrent/{ASSEMBLY}.TOPUP.recurrent.tsv", ASSEMBLY=ASSEMBLIES),
        expand("csv/orthogonal_support/all/{SAMPLE}.{ASSEMBLY}.{TECH}.TOPUP.tsv", SAMPLE=ALL_SAMPLES, ASSEMBLY=ASSEMBLIES, TECH=["hifi"])
