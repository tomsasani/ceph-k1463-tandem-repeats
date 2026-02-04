import pandas as pd
import glob

wildcard_constraints:
    SAMPLE = r"[0-9]{4,6}",
    ASSEMBLY = "GRCh38|CHM13v2",
    TECH = "element|hifi|illumina",
    TOPUP = "TOPUP|ORIGINAL",
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
    "GRCh38": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.bed.gz",
    "CHM13v2": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/chm13v2.0_maskedY_rCRS.palladium-v1.0.trgt.bed.gz",
    # "GRCh38": "data/gangstr.trgt_formatted.bed.gz"
    }

HPRC_PREF = f"/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/"

HPRC_SAMPLES = []
for fh in glob.glob(HPRC_PREF + "GRCh38_v1.0_50bp_merge/1.1.2-69937d83/hprc/*.vcf.gz"):
    sample_id = fh.split("/")[-1].split("_")[0]
    if sample_id == "hprc": continue
    HPRC_SAMPLES.append(sample_id)


def get_complete_bams(wildcards):
    """
    return the path to the 'complete' BAM file for a given sample.
    if the sample is one of the four with top-up sequencing, return a path
    that requires us to re-align and merge the updated top-up data. otherwise
    return the original path to the HiFi BAM
    """
    if wildcards.USE_NEW_BAM == "TOPUP" and wildcards.SAMPLE in ("200100", "2189", "2216", "200080"):
        return f"{NEW_PATH}/merged/{wildcards.ASSEMBLY}/{wildcards.SAMPLE}.merged.bam"
    else:
        assembly_adj = wildcards.ASSEMBLY.split('v2')[0]
        sample_id = SMP2ALT[wildcards.SAMPLE]
        return "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams/{0}/{1}.{2}.haplotagged.bam".format(assembly_adj, sample_id, assembly_adj.lower() if 'CHM' in wildcards.ASSEMBLY else assembly_adj,)


def get_bam_for_validation(sample, assembly, tech, use_new_bam):

    assembly_adj = assembly.split('v2')[0]
    sample_id = SMP2ALT[sample] if tech in ("ont", "hifi") else sample

    tech2path = {
        "ont": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/ont-bams/{0}/{1}.minimap2.bam".format(assembly_adj, sample_id,), 
        "element": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/element/{0}_bams/{1}-E.{0}.merged.sort.bam".format(assembly_adj, sample_id,),
        "illumina": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/CEPH/cram/{0}.cram".format(sample_id),
        "hifi": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams/{0}/{1}.{2}.haplotagged.bam".format(assembly_adj, sample_id, assembly_adj.lower() if 'CHM' in assembly else assembly_adj,),
        }

    if sample in ("2216", "2189", "200080", "200100") and use_new_bam == "TOPUP":
        tech2path.update({"hifi": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/UW_PB_HiFi/topoff_2025/merged/{0}/{1}.merged.bam".format(assembly_adj, sample)})
    
    return tech2path[tech]


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
    return ASSEMBLY2CATALOG[wildcards.ASSEMBLY]


# create global dictionaries we'll use
PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str})
ped["sex"] = ped["suffix"].apply(lambda s: "male" if ("S" in s or s == "F") else "female" )

SMP2SEX = dict(zip(ped["sample_id"], ped["sex"]))
SMP2DAD = dict(zip(ped["sample_id"], ped["paternal_id"]))
SMP2MOM = dict(zip(ped["sample_id"], ped["maternal_id"]))
SMP2ALT = dict(zip(ped["sample_id"], ped["alt_sample_id"]))
SMP2SUFF = dict(zip(ped["sample_id"], ped["suffix"]))


CHILDREN = ped[(ped["paternal_id"] != "missing")]["sample_id"].to_list()
ALL_SAMPLES = ped["sample_id"].to_list()
ASSEMBLIES = ["CHM13v2"]


include: "rules/alignment.smk"
include: "rules/snvs.smk"
include: "rules/trgt.smk"
include: "rules/phasing.smk"
include: "rules/trgt_denovo.smk"
include: "rules/denovo_mutation_filtering.smk"


rule all:
    input:
        expand("csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.TOPUP.tsv", SAMPLE=CHILDREN, ASSEMBLY=ASSEMBLIES),
        expand("csv/denominators/{SAMPLE}.{ASSEMBLY}.TOPUP.denominator.tsv", SAMPLE=CHILDREN, ASSEMBLY=ASSEMBLIES),
        expand("csv/annotated/{SAMPLE}.{ASSEMBLY}.{TECH}.TOPUP.tsv", SAMPLE=CHILDREN, ASSEMBLY=ASSEMBLIES, TECH=["element", "hifi"]),
        expand("csv/recurrent/{ASSEMBLY}.TOPUP.recurrent.tsv", ASSEMBLY=ASSEMBLIES),
        expand("csv/orthogonal_support/all/{SAMPLE}.{ASSEMBLY}.{TECH}.TOPUP.tsv", SAMPLE=ALL_SAMPLES, ASSEMBLY=ASSEMBLIES, TECH=["hifi"]),
        expand("csv/phasing/{SAMPLE}.{ASSEMBLY}.TOPUP.phased.3gen.tsv", SAMPLE=["2216", "2189"], ASSEMBLY=ASSEMBLIES),
        # "csv/hprc/combined.CHM13v2.heterozygosity.tsv"