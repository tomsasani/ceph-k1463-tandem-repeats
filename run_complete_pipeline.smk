include: "rules/alignment.smk"
include: "rules/snvs.smk"
include: "rules/trgt.smk"
include: "rules/phasing.smk"
include: "rules/trgt_denovo.smk"
include: "rules/denovo_mutation_filtering.smk"


rule all:
    input:
        "csv/filtered_and_merged/200101.CHM13v2.TOPUP.tsv",
        #"csv/raw_denovos/2189.CHM13v2.TOPUP.chr22.trgt-denovo.csv"