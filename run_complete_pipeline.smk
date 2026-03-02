import pandas as pd
import glob

wildcard_constraints:
    SAMPLE = r"(HG)?(NA)?[0-9]{4,6}",
    ASSEMBLY = "GRCh38|CHM13v2",
    TECH = "element|hifi|illumina",
    USE_NEW_BAM = "TOPUP|ORIGINAL",
    CHROM = r"chr[0-9]{1,2}|chrX|chrY"


CHROMS = list(map(str, range(1, 23)))
CHROMS = [f"chr{c}" for c in CHROMS]
CHROMS.extend(["chrX", "chrY"])

NEW_PATH = "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/UW_PB_HiFi/topoff_2025"
CUR_PREF = "/scratch/ucgd/lustre-labs/quinlan/u1006375/CEPH-K1463-TandemRepeats/"

ASSEMBLY2REF = {
    "GRCh38": "/scratch/ucgd/lustre/common/data/Reference/homo_sapiens/GRCh38/primary_assembly_decoy_phix.fa",
    "CHM13v2": "/scratch/ucgd/lustre/common/data/Reference/homo_sapiens/CHM13v2.0/primary_assembly_decoy_phix.fa"
}

ASSEMBLY2CATALOG = {
    # "GRCh38": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.bed.gz",
    "CHM13v2": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/chm13v2.0_maskedY_rCRS.palladium-v1.0.trgt.bed.gz",
    "GRCh38": "data/gangstr.trgt_formatted.bed.gz"
}

ASSEMBLY2ANNOTATIONS = {
    # "GRCh38": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.annotations.bed.gz",
    "GRCh38": "data/gangstr.annotations.bed.gz",
    "CHM13v2": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/chm13v2.0_maskedY_rCRS.palladium-v1.0.trgt.annotations.bed.gz",
}

HPRC_PREF = f"/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/GRCh38_v1.0_50bp_merge/1.1.2-69937d83/hprc"

HPRC_SAMPLES = []
for fh in glob.glob(f"{HPRC_PREF}/*.vcf.gz"):
    sample_id = fh.split("/")[-1].split("_")[0]
    if sample_id == "hprc": continue
    HPRC_SAMPLES.append(sample_id)

BAM_MAP = pd.read_csv("K1463_bam_mapping.tsv", sep="\t", dtype={"sample_id": str})

TOPPED_UP_SAMPLES = ["2216", "2189", "200100", "200080"]

def get_complete_bams(wildcards):
    """
    return the path to the 'complete' BAM file for a given sample.
    if the sample is one of the four with top-up sequencing, return a path
    that requires us to re-align and merge the updated top-up data. otherwise
    return the original path to the HiFi BAM
    """
    try:
        technology = wildcards.TECH
    except AttributeError:
        technology = "hifi"

    path_key = "topup_path" if wildcards.USE_NEW_BAM == "TOPUP" and wildcards.SAMPLE in TOPPED_UP_SAMPLES and technology == "hifi" else "path"

    return BAM_MAP[(BAM_MAP["sample_id"] == wildcards.SAMPLE) & \
            (BAM_MAP["assembly"] == wildcards.ASSEMBLY) & \
            (BAM_MAP["tech"] == technology)][path_key].to_list()[0]

def get_bam_for_validation(sample, assembly, tech, use_new_bam):

    path_key = "topup_path" if use_new_bam == "TOPUP" and sample in TOPPED_UP_SAMPLES and tech == "hifi" else "path"
    
    return BAM_MAP[(BAM_MAP["sample_id"] == sample) & \
            (BAM_MAP["assembly"] == assembly) & \
            (BAM_MAP["tech"] == tech)][path_key].to_list()[0]


def get_children_vcfs(wildcards):
    if wildcards.SAMPLE == "2216":
        children = ["200081", "200082", "200084", "200085", "200086", "200087"]
    elif wildcards.SAMPLE == "2189":
        children = ["200101", "200102", "200103", "200104", "200105", "200106"]
    else:
        children = []
    
    return [f"data/trgt/phased/{s}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.phased.vcf.gz" for s in children]


def get_grandparent_vcfs(wildcards):
    if wildcards.SAMPLE.startswith("200"):
        grandparents = ["2209", "2188"]
    else:
        grandparents = ["2281", "2280", "2213", "2214"]
        
    return [f"data/trgt/phased/{s}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.phased.vcf.gz" for s in grandparents]


def get_annotation_fh(wildcards):
    return ASSEMBLY2ANNOTATIONS[wildcards.ASSEMBLY]


# create global dictionaries we'll use
PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str})
ped["sex"] = ped["suffix"].apply(lambda s: "male" if ("S" in s or s == "F") else "female" )

SMP2SEX = dict(zip(ped["sample_id"], ped["sex"]))
SMP2DAD = dict(zip(ped["sample_id"], ped["paternal_id"]))
SMP2MOM = dict(zip(ped["sample_id"], ped["maternal_id"]))
SMP2ALT = dict(zip(ped["sample_id"], ped["alt_sample_id"]))
SMP2SUFF = dict(zip(ped["sample_id"], ped["suffix"]))


CHILDREN = ped[(ped["paternal_id"] != "missing") & (~ped["paternal_id"].isin(["2281", "2214"]))]["sample_id"].to_list()
ALL_SAMPLES = ped["sample_id"].to_list()
ASSEMBLIES = ["GRCh38"]

TRIDS = ["chr8_2623352_2623487_trsolve"]#, "chr1_54510591_54510710_trsolve"]

USE_NEW_BAMS = ["TOPUP"]

include: "rules/alignment.smk"
include: "rules/snvs.smk"
include: "rules/trgt.smk"
include: "rules/phasing.smk"
include: "rules/trgt_denovo.smk"
include: "rules/denovo_mutation_filtering.smk"
include: "rules/decomposition.smk"


rule all:
    input:
        expand("csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.tsv", SAMPLE=CHILDREN, ASSEMBLY=ASSEMBLIES, USE_NEW_BAM=USE_NEW_BAMS),
        expand("csv/denominators/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.denominator.tsv", SAMPLE=CHILDREN, ASSEMBLY=ASSEMBLIES, USE_NEW_BAM=USE_NEW_BAMS),
        expand("csv/annotated/{SAMPLE}.{ASSEMBLY}.{TECH}.{USE_NEW_BAM}.tsv", SAMPLE=CHILDREN, ASSEMBLY=ASSEMBLIES, TECH=["element", "hifi"], USE_NEW_BAM=USE_NEW_BAMS),
        expand("csv/recurrent/{ASSEMBLY}.{USE_NEW_BAM}.recurrent.tsv", ASSEMBLY=ASSEMBLIES, USE_NEW_BAM=USE_NEW_BAMS),
        expand("csv/orthogonal_support/all/{SAMPLE}.{ASSEMBLY}.{TECH}.{USE_NEW_BAM}.tsv", SAMPLE=ALL_SAMPLES, ASSEMBLY=ASSEMBLIES, TECH=["hifi"], USE_NEW_BAM=USE_NEW_BAMS),
        expand("csv/phasing/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.3gen.tsv", SAMPLE=["2216", "2189"], ASSEMBLY=ASSEMBLIES, USE_NEW_BAM=USE_NEW_BAMS),
        # expand("csv/hprc/combined.{ASSEMBLY}.{USE_NEW_BAM}.heterozygosity.tsv", ASSEMBLY=ASSEMBLIES, USE_NEW_BAM=USE_NEW_BAMS),
        # expand("trviz/{COHORT}/{TRID}.GRCh38.key.tsv", TRID=TRIDS, COHORT=["ceph", "hprc"])