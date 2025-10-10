import pandas as pd


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


# create global dictionaries we'll use
PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str})

ped["sex"] = ped["suffix"].apply(lambda s: "male" if ("S" in s or s == "F") else "female" )
SMP2SEX = dict(zip(ped["sample_id"], ped["sex"]))


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


rule convert_gangstr_to_trgt:
    input:
        gangstr = "data/hg38_ver17.bed.gz"
    output:
        trgt = "data/gangstr.trgt_formatted.bed"
    script:
        "scripts/convert_gangstr_to_trgt.py"


rule compress_trgt:
    input:
        "data/gangstr.trgt_formatted.bed"
    output:
        "data/gangstr.trgt_formatted.bed.gz"
    shell:
        """
        gzip {input}
        """


rule create_chrom_bed:
    input:
        repeats = lambda wildcards: ASSEMBLY2CATALOG[wildcards.ASSEMBLY]
    output:
        "data/catalogs/{CHROM}.{ASSEMBLY}.catalog.bed"
    shell:
        """
        zgrep -w {wildcards.CHROM} {input.repeats} > {output}
        """


rule create_chrom_ref:
    input:
        reference = lambda wildcards: ASSEMBLY2REF[wildcards.ASSEMBLY],
    output: contig = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz"
    shell:
        """
        module load samtools
        
        samtools faidx {input.reference} {wildcards.CHROM} | bgzip > {output.contig}
        """


rule index_chrom_ref:
    input:
        contig = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz"
    output: idx = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz.fai"
    shell:
        """
        module load samtools
        
        samtools faidx {input.contig}
        """


rule run_trgt:
    input:
        bam = lambda wildcards: get_complete_bams(wildcards),
        bam_idx = lambda wildcards: get_complete_bams(wildcards) + ".bai",
        reference = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz",
        trgt_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/src/trgt-v4.0.0-x86_64-unknown-linux-gnu/trgt",
        repeats = "data/catalogs/{CHROM}.{ASSEMBLY}.catalog.bed"
    output:
        vcf = temp("data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.vcf.gz"),
        bam = temp("data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.spanning.bam")
    params:
        karyotype_cmd = lambda wildcards: "--karyotype XY" if SMP2SEX[wildcards.SAMPLE] == "male" and wildcards.CHROM in ("chrX", "chrY") else "",
        output_prefix = lambda wildcards: f"data/trgt/per-chrom/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}"
    threads: 16
    resources:
        mem_mb = 64_000
    script:
        "bash_scripts/run_trgt.sh"


rule sort_chrom_vcfs:
    input: 
        vcf = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.vcf.gz",
    output:
        vcf = CUR_PREF + "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz",
    script:
        "bash_scripts/sort_trgt_vcf.sh"


rule index_trgt_vcfs:
    input:
        vcf = CUR_PREF + "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz"
    output:
        CUR_PREF + "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz.tbi"
    script:
        "bash_scripts/index_vcf.sh"


rule sort_bam:
    input: bam = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.spanning.bam"
    output: bam = CUR_PREF + "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.spanning.sorted.bam"
    script:
        "bash_scripts/sort_bam.sh"


rule index_bam:
    input: bam = CUR_PREF + "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.spanning.sorted.bam"
    output: CUR_PREF + "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.spanning.sorted.bam.bai"
    script:
        "bash_scripts/index_bam.sh"