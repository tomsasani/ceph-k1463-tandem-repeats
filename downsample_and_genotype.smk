import pandas as pd 
from typing import Dict


CHROMS = [f"chr{c}" for c in range(1, 23)] + ["chrX"]

PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str,})
ped["sex"] = ped["suffix"].apply(lambda s: "male" if ("S" in s or s == "F") else "female" )
SMP2SEX = dict(zip(ped["sample_id"], ped["sex"]))
SAMPLES = ped[ped["paternal_id"] == "2209"]["sample_id"].to_list()
SMP2ALT = dict(zip(ped["sample_id"], ped["alt_sample_id"]))
SMP2SUFF = dict(zip(ped["sample_id"].to_list(), ped["suffix"].to_list()))
SMP2DAD = dict(zip(ped["sample_id"].to_list(), ped["paternal_id"].to_list()))
SMP2MOM = dict(zip(ped["sample_id"].to_list(), ped["maternal_id"].to_list()))

SAMPLES = ped[ped["paternal_id"] == "2209"]["sample_id"].to_list()
PARENTS = ["2209", "2188"]
CHILDREN = ped[ped["paternal_id"] != "missing"]["sample_id"].to_list()
ALL_SAMPLES = ped["sample_id"].to_list()

CUR_PREF = "/scratch/ucgd/lustre-labs/quinlan/u1006375/CEPH-K1463-TandemRepeats/"


ASSEMBLY2CATALOG = {"GRCh38": "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.bed.gz",
                   "CHM13v2": "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/chm13v2.0_maskedY_rCRS.palladium-v1.0.trgt.bed.gz"}

ASSEMBLY2REF = {"GRCh38": "/scratch/ucgd/lustre/common/data/Reference/homo_sapiens/GRCh38/primary_assembly_decoy_phix.fa", 
                "CHM13v2": "/scratch/ucgd/lustre/common/data/Reference/homo_sapiens/CHM13v2.0/primary_assembly_decoy_phix.fa"}


def get_complete_bams(wildcards):
    """
    return the path to the 'complete' BAM file for a given sample.
    if the sample is one of the four with top-up sequencing, return a path
    that requires us to re-align and merge the updated top-up data. otherwise
    return the original path to the HiFi BAM
    """

    # define the sample we're getting the BAM for
    try:
        if wildcards.MEMBER == "kid":
            sample = wildcards.SAMPLE
        elif wildcards.MEMBER == "mom":
            sample = SMP2MOM[wildcards.SAMPLE]
        elif wildcards.MEMBER == "dad":
            sample = SMP2DAD[wildcards.SAMPLE]
    except AttributeError:
        sample = wildcards.SAMPLE

    return BAM_MAP[(BAM_MAP["sample_id"] == sample) & \
            (BAM_MAP["assembly"] == wildcards.ASSEMBLY) & \
            (BAM_MAP["tech"] == "hifi")]["path"].to_list()[0]


def get_annotation_fh(wildcards):
    if wildcards.ASSEMBLY == "GRCh38":
        return "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.annotations.bed.gz"
    elif wildcards.ASSEMBLY == "CHM13v2":
        return "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/chm13v2.0_maskedY_rCRS.palladium-v1.0.trgt.annotations.bed.gz"


def get_grandparent_vcfs(wildcards):
    if wildcards.SAMPLE.startswith("200"):
        grandparents = ["2209", "2188"]
    else:
        grandparents = ["2281", "2280", "2213", "2214"]
        
    return [f"data/trgt/phased/{s}.{wildcards.ASSEMBLY}.TOPUP.phased.vcf.gz" for s in grandparents]


wildcard_constraints:
    SAMPLE = "[0-9]+",
    ASSEMBLY = "CHM13v2|GRCh38",
    CHROM = "chr[0-9]{1,2}|chrX|chrY",


BAM_MAP = pd.read_csv("K1463_bam_mapping.tsv", sep="\t", dtype={"sample_id": str})


KID_TARGETS = [50, 50, 50, 50, 50]
DAD_TARGETS = [10, 30, 50, 50, 50]
MOM_TARGETS = [50, 50, 50, 30, 10]

rule all:
    input:
        expand("downsampling/csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.tsv", zip, 
            SAMPLE = ["2217"] * len(KID_TARGETS), 
            ASSEMBLY = ["CHM13v2"] * len(KID_TARGETS), 
            KID_TARGET=KID_TARGETS,
            MOM_TARGET=MOM_TARGETS, 
            DAD_TARGET=DAD_TARGETS,
        )


rule run_mosdepth:
    input:
        mosdepth = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/mosdepth",
        bam = get_complete_bams,
    output:
        summary = "downsampling/data/mosdepth/{SAMPLE}.{ASSEMBLY}.mosdepth.summary.txt"
    shell:
        """
        {input.mosdepth} --by 1000 --fast-mode -n downsampling/data/mosdepth/{wildcards.SAMPLE}.{wildcards.ASSEMBLY} {input.bam}
        """


rule output_mosdepth_summary:
    input:
        kid_depth = lambda wildcards: f"downsampling/data/mosdepth/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.mosdepth.summary.txt",
        mom_depth = lambda wildcards: f"downsampling/data/mosdepth/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.mosdepth.summary.txt",
        dad_depth = lambda wildcards: f"downsampling/data/mosdepth/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.mosdepth.summary.txt"
    output:
        mean_depths = "downsampling/csv/{SAMPLE}.{ASSEMBLY}.depth.txt",
    run:
        import pandas as pd

        kid_depth = pd.read_csv(input.kid_depth, sep="\t")
        mean_kid_depth = kid_depth[kid_depth["chrom"] == "total"]["mean"].values[0]

        dad_depth = pd.read_csv(input.dad_depth, sep="\t")
        mean_dad_depth = dad_depth[dad_depth["chrom"] == "total"]["mean"].values[0]

        mom_depth = pd.read_csv(input.mom_depth, sep="\t")
        mean_mom_depth = mom_depth[mom_depth["chrom"] == "total"]["mean"].values[0]
        
        with open(output.mean_depths, "w") as outfh:
            print ("\t".join(["dad_dp_mean", "mom_dp_mean", "kid_dp_mean"]), file=outfh)
            print ("\t".join(list(map(str, [mean_dad_depth, mean_mom_depth, mean_kid_depth]))), file=outfh)
        outfh.close()


rule downsample_bam:
    input:
        py_script = "rules/scripts/create_downsampling_cmd.py",
        mean_depths = "downsampling/csv/{SAMPLE}.{ASSEMBLY}.depth.txt",
        bam = get_complete_bams
    output: 
        bam = "downsampling/data/bam/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.subsampled.bam",
    threads: 16
    shell:
        """
        module load sambamba

        echo {wildcards.MEMBER}

        CMD=$(python {input.py_script} --depth {input.mean_depths} \
                                 --input_bam {input.bam} \
                                 --output_bam {output.bam} \
                                 -target_depth {wildcards.TARGET} \
                                 -nthreads 16 \
                                 -member {wildcards.MEMBER})

        echo $CMD
        echo $CMD | sh
        """


rule index_downsampled_bam:
    input:
        bam = "downsampling/data/bam/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.subsampled.bam"
    output:
        bam_idx = "downsampling/data/bam/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.subsampled.bam.bai"
    script:
        "rules/bash_scripts/index_bam.sh"


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
        bam = "downsampling/data/bam/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.subsampled.bam",
        bam_idx = "downsampling/data/bam/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.subsampled.bam.bai",
        reference = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz",
        trgt_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/src/trgt-v4.0.0-x86_64-unknown-linux-gnu/trgt",
        repeats = "data/catalogs/{CHROM}.{ASSEMBLY}.catalog.bed"
    output:
        vcf = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.vcf.gz",
        bam = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.spanning.bam"
    params:
        karyotype_cmd = lambda wildcards: "--karyotype XY" if SMP2SEX[wildcards.SAMPLE] == "male" and wildcards.CHROM in ("chrX", "chrY") else "",
        output_prefix = lambda wildcards: f"downsampling/data/trgt/per-chrom/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.{wildcards.MEMBER}.{wildcards.TARGET}.{wildcards.CHROM}"
    threads: 16
    resources:
        mem_mb = 64_000
    script:
        "rules/bash_scripts/run_trgt.sh"


rule sort_chrom_vcfs:
    input: 
        vcf = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.vcf.gz",
    output:
        vcf = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.sorted.vcf.gz",
    script:
        "rules/bash_scripts/sort_trgt_vcf.sh"

rule index_chrom_vcfs:
    input: 
        vcf = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.sorted.vcf.gz",
    output: "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.sorted.vcf.gz.tbi",
    script:
        "rules/bash_scripts/index_vcf.sh"


rule sort_bam:
    input: bam = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.spanning.bam"
    output: bam = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.spanning.sorted.bam"
    script:
        "rules/bash_scripts/sort_bam.sh"


rule index_bam:
    input: bam = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.spanning.sorted.bam"
    output: "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.spanning.sorted.bam.bai"
    script:
        "rules/bash_scripts/index_bam.sh"


def get_alt_id(wildcards):
    if wildcards.MEMBER == "kid":
        return SMP2ALT[wildcards.SAMPLE]
    elif wildcards.MEMBER == "mom":
        return SMP2ALT[SMP2MOM[wildcards.SAMPLE]]
    elif wildcards.MEMBER == "dad":
        return SMP2ALT[SMP2DAD[wildcards.SAMPLE]]

rule call_snvs:
    input:
        ref = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz",
        ref_idx = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz.fai",
        bam = "downsampling/data/bam/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.subsampled.bam",
        bam_idx = "downsampling/data/bam/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.subsampled.bam.bai",
        sif = "deepvariant_1.9.0.sif"
    output:
        vcf = "downsampling/data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.vcf.gz",
        gvcf = "downsampling/data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.g.vcf.gz",
    threads: 8
    params:
        alt_sample_id = lambda wildcards: get_alt_id(wildcards)
    resources:
        mem_mb = 64_000,
    script:
        "rules/bash_scripts/run_deepvariant.sh"


rule sort_snv_chrom_vcfs:
    input: 
        vcf = "downsampling/data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.vcf.gz",
    output: vcf = "downsampling/data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.sorted.vcf.gz"
    script:
       "rules/bash_scripts/sort_snv_vcf.sh"


rule index_snv_chrom_vcfs:
    input: vcf = "downsampling/data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.sorted.vcf.gz"
    output: "downsampling/data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.sorted.vcf.gz.tbi"
    threads: 4
    script:
        "rules/bash_scripts/index_vcf.sh"


rule combine_trio_chrom_vcfs:
    input:
        sif = "glnexus_v1.2.7.sif",
        kid_snp_gvcf = "downsampling/data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.kid.{KID_TARGET}.{CHROM}.g.vcf.gz",
        mom_snp_gvcf = "downsampling/data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.mom.{MOM_TARGET}.{CHROM}.g.vcf.gz",
        dad_snp_gvcf = "downsampling/data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.dad.{DAD_TARGET}.{CHROM}.g.vcf.gz",
    output: "downsampling/data/vcf/trios/per-chrom/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.{CHROM}.trio.vcf.gz"
    params:
        gl_nexus_prefix = lambda wildcards: f"gl_nexus_dbs/{wildcards.SAMPLE}_{wildcards.ASSEMBLY}_{wildcards.KID_TARGET}_{wildcards.MOM_TARGET}_{wildcards.DAD_TARGET}_{wildcards.CHROM}"
    threads: 4
    resources:
        mem_mb = 32_000
    script:
        "rules/bash_scripts/combine_trio_vcf.sh"


rule merge_trio_vcfs:
    input:
        vcfs = expand("downsampling/data/vcf/trios/per-chrom/{{SAMPLE}}.{{ASSEMBLY}}.{{KID_TARGET}}.{{MOM_TARGET}}.{{DAD_TARGET}}.{CHROM}.trio.vcf.gz", CHROM=CHROMS)
    output: "downsampling/data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.vcf.gz"
    threads: 8
    shell:
        """
        module load bcftools
        
        bcftools concat {input.vcfs} | bcftools view --threads {threads} | bgzip > {output}
        """


rule index_merged_vcf:
    input:
        vcf = "downsampling/data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.vcf.gz"
    output: "downsampling/data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.vcf.gz.tbi"
    script:
        "rules/bash_scripts/index_vcf.sh"


rule run_hiphase:
    input:
        snv_vcf = "downsampling/data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.sorted.vcf.gz",
        snv_vcf_idx = "downsampling/data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.sorted.vcf.gz.tbi",
        str_vcf = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.sorted.vcf.gz",
        str_vcf_idx = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.sorted.vcf.gz.tbi",
        bam = "downsampling/data/bam/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.subsampled.bam",
        reference = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz"
    output:
        snv_vcf = "downsampling/data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.phased.vcf.gz",
        str_vcf = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.{CHROM}.phased.vcf.gz",
    threads: 4
    script:
        "rules/bash_scripts/run_hiphase.sh"


rule combine_phased_snv_vcfs:
    input:
        vcfs = expand("downsampling/data/vcf/snv/per-chrom/{{SAMPLE}}.{{ASSEMBLY}}.{{MEMBER}}.{{TARGET}}.{CHROM}.phased.vcf.gz", CHROM=CHROMS)
    output: "downsampling/data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.phased.vcf.gz"
    threads: 8
    shell:
        """
        module load bcftools
        
        bcftools concat {input.vcfs} | bcftools view --threads {threads} -f 'PASS' | bgzip > {output}
        """


rule index_phased_snv_vcfs:
    input: vcf = "downsampling/data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.phased.vcf.gz"
    output: "downsampling/data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.phased.vcf.gz.tbi"
    script:
        "rules/bash_scripts/index_vcf.sh"


rule combine_phased_trgt_vcfs:
    input:
        vcfs = expand("downsampling/data/trgt/per-chrom/{{SAMPLE}}.{{ASSEMBLY}}.{{MEMBER}}.{{TARGET}}.{CHROM}.phased.vcf.gz", CHROM=CHROMS)
    output: "downsampling/data/trgt/phased/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.phased.vcf.gz"
    threads: 8
    shell:
        """
        module load bcftools
        
        bcftools concat {input.vcfs} | bcftools view --threads {threads} | bgzip > {output}
        """


rule index_phased_trgt_vcfs:
    input: vcf = "downsampling/data/trgt/phased/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.phased.vcf.gz"
    output: "downsampling/data/trgt/phased/{SAMPLE}.{ASSEMBLY}.{MEMBER}.{TARGET}.phased.vcf.gz.tbi"
    script:
        "rules/bash_scripts/index_vcf.sh"


rule run_trgt_denovo:
    input:
        reference = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz",
        repeats = "data/catalogs/{CHROM}.{ASSEMBLY}.catalog.bed",
        trgt_denovo_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/trgt-denovo",

        kid_vcf = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.kid.{KID_TARGET}.{CHROM}.sorted.vcf.gz",
        kid_vcf_idx = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.kid.{KID_TARGET}.{CHROM}.sorted.vcf.gz.tbi",

        mom_vcf = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.mom.{MOM_TARGET}.{CHROM}.sorted.vcf.gz",
        mom_vcf_idx = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.mom.{MOM_TARGET}.{CHROM}.sorted.vcf.gz.tbi",

        dad_vcf = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.dad.{DAD_TARGET}.{CHROM}.sorted.vcf.gz",
        dad_vcf_idx = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.dad.{DAD_TARGET}.{CHROM}.sorted.vcf.gz.tbi",

        kid_bam = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.kid.{KID_TARGET}.{CHROM}.spanning.sorted.bam",
        kid_bam_idx = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.kid.{KID_TARGET}.{CHROM}.spanning.sorted.bam.bai",

        mom_bam = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.mom.{MOM_TARGET}.{CHROM}.spanning.sorted.bam",
        mom_bam_idx = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.mom.{MOM_TARGET}.{CHROM}.spanning.sorted.bam.bai",

        dad_bam = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.dad.{DAD_TARGET}.{CHROM}.spanning.sorted.bam",
        dad_bam_idx = "downsampling/data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.dad.{DAD_TARGET}.{CHROM}.spanning.sorted.bam.bai",
    output:
        output = "downsampling/csv/denovos/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.{CHROM}.denovo.csv"
    threads: 32
    resources:
        mem_mb = 32_000,
        runtime = 720,
        slurm_account = "quinlan-rw",
        slurm_partition = "quinlan-rw"
    params:
        output_dir = lambda wildcards: f"downsampling_tmpdir/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.{wildcards.KID_TARGET}.{wildcards.MOM_TARGET}.{wildcards.DAD_TARGET}.{wildcards.CHROM}",
        kid_pref = lambda wildcards: f"downsampling/data/trgt/per-chrom/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.kid.{wildcards.KID_TARGET}.{wildcards.CHROM}",
        mom_pref = lambda wildcards: f"downsampling/data/trgt/per-chrom/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.mom.{wildcards.MOM_TARGET}.{wildcards.CHROM}",
        dad_pref = lambda wildcards: f"downsampling/data/trgt/per-chrom/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.dad.{wildcards.DAD_TARGET}.{wildcards.CHROM}",
        global_pref =  CUR_PREF,

    script:
        "rules/bash_scripts/run_trgt_denovo.sh"


rule combine_trgt_denovo:
    input:
        fhs = expand("downsampling/csv/denovos/{{SAMPLE}}.{{ASSEMBLY}}.{{KID_TARGET}}.{{MOM_TARGET}}.{{DAD_TARGET}}.{CHROM}.denovo.csv", CHROM=CHROMS)
    output: fh = "downsampling/csv/combined/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.denovo.csv"
    shell:
        """
        
        grep '^chrom' {input.fhs[0]} > {output.fh}

        cat {input.fhs} | grep -v '^chrom'  >> {output.fh}
        """


rule prefilter_denovos:
    input:
        ped = "tr_validation/data/file_mapping.csv",
        kid_mutation_df = "downsampling/csv/combined/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.denovo.csv",
        annotations = lambda wildcards: get_annotation_fh(wildcards),
        utils = "rules/scripts/utils.py"
    output: 
        fh = "downsampling/csv/prefiltered/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.denovo.csv"
    params:
        filtering_mode = "strict"
    script:
        "rules/scripts/filter_mutation_dataframe.py"


rule annotate_with_informative_sites:
    input:
        cohort_snp_vcf = "downsampling/data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.vcf.gz",
        cohort_snp_vcf_idx = "downsampling/data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.vcf.gz.tbi",
        kid_phased_snp_vcf = "downsampling/data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.kid.{KID_TARGET}.phased.vcf.gz",
        kid_phased_snp_vcf_idx = "downsampling/data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.kid.{KID_TARGET}.phased.vcf.gz.tbi",
        kid_mutation_df = "downsampling/csv/prefiltered/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.denovo.csv",
        kid_phased_str_vcf = "downsampling/data/trgt/phased/{SAMPLE}.{ASSEMBLY}.kid.{KID_TARGET}.phased.vcf.gz",
        kid_phased_str_vcf_idx = "downsampling/data/trgt/phased/{SAMPLE}.{ASSEMBLY}.kid.{KID_TARGET}.phased.vcf.gz.tbi",
    output:
        out = "downsampling/csv/phasing/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.annotated.2gen.tsv"
    params:
        dad_id = lambda wildcards: SMP2ALT[SMP2DAD[wildcards.SAMPLE]],
        mom_id = lambda wildcards: SMP2ALT[SMP2MOM[wildcards.SAMPLE]],
        focal_alt_id = lambda wildcards: SMP2ALT[wildcards.SAMPLE],
        slop = 500_000
    script:
        "rules/scripts/annotate_with_informative_sites.py"


rule annotate_with_parental_haplotype:
    input:
        annotated_dnms = "downsampling/csv/phasing/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.annotated.2gen.tsv",
        dad_phased_str_vcf = lambda wildcards: f"downsampling/data/trgt/phased/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.dad.{wildcards.DAD_TARGET}.phased.vcf.gz",
        dad_phased_str_vcf_idx = lambda wildcards: f"downsampling/data/trgt/phased/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.dad.{wildcards.DAD_TARGET}.phased.vcf.gz.tbi",
        dad_phased_snv_vcf = lambda wildcards: f"downsampling/data/vcf/snv/phased/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.dad.{wildcards.DAD_TARGET}.phased.vcf.gz",
        dad_phased_snv_vcf_idx = lambda wildcards: f"downsampling/data/vcf/snv/phased/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.dad.{wildcards.DAD_TARGET}.phased.vcf.gz.tbi",
        mom_phased_str_vcf = lambda wildcards: f"downsampling/data/trgt/phased/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.mom.{wildcards.MOM_TARGET}.phased.vcf.gz",
        mom_phased_str_vcf_idx = lambda wildcards: f"downsampling/data/trgt/phased/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.mom.{wildcards.MOM_TARGET}.phased.vcf.gz.tbi",
        mom_phased_snv_vcf = lambda wildcards: f"downsampling/data/vcf/snv/phased/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.mom.{wildcards.MOM_TARGET}.phased.vcf.gz",
        mom_phased_snv_vcf_idx = lambda wildcards: f"downsampling/data/vcf/snv/phased/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.mom.{wildcards.MOM_TARGET}.phased.vcf.gz.tbi",
    output:
        out = "downsampling/csv/phasing/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.haplotyped.2gen.tsv"
    params:
        dad_id = lambda wildcards: SMP2ALT[SMP2DAD[wildcards.SAMPLE]],
        mom_id = lambda wildcards: SMP2ALT[SMP2MOM[wildcards.SAMPLE]],
    script:
        "rules/scripts/annotate_with_parental_haplotype.py"


rule phase:
    input:
        annotated_dnms = "downsampling/csv/phasing/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.haplotyped.2gen.tsv",
        py_script = "rules/scripts/phase_by_ps.py",
    output: "downsampling/csv/phasing/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.phased.2gen.tsv"
    shell:
        """
        python {input.py_script} --annotated_dnms {input.annotated_dnms} \
                                 --out {output}
        """


rule annotate_with_parental_alleles:
    input:
        annotated_dnms = "downsampling/csv/phasing/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.phased.2gen.tsv",
        dad_phased_str_vcf = "downsampling/data/trgt/phased/{SAMPLE}.{ASSEMBLY}.dad.{DAD_TARGET}.phased.vcf.gz",
        mom_phased_str_vcf = "downsampling/data/trgt/phased/{SAMPLE}.{ASSEMBLY}.mom.{MOM_TARGET}.phased.vcf.gz",
        kid_phased_str_vcf = "downsampling/data/trgt/phased/{SAMPLE}.{ASSEMBLY}.kid.{KID_TARGET}.phased.vcf.gz",
    output:
        out = "downsampling/csv/phasing/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.allele_sequences.tsv"
    script:
        "rules/scripts/annotate_with_parental_allele_sequences.py"


rule add_transmission_evidence:
    input:
        mutations = "downsampling/csv/prefiltered/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.denovo.csv",
        other_vcfs = lambda wildcards: []#get_children_vcfs
    output: fh = "downsampling/csv/transmission/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.tsv"
    params:
        generation_to_query = "children"
    script:
        "rules/scripts/assess_presence_in_cohort.py"


rule add_grandparental_evidence:
    input:
        mutations = "downsampling/csv/prefiltered/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.denovo.csv",
        other_vcfs = get_grandparent_vcfs
    output: fh = "downsampling/csv/grandparents/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.tsv"
    params:
        generation_to_query = "grandparents"
    script:
        "rules/scripts/assess_presence_in_cohort.py"


# take all of the de novo mutation metadata and merge into one combined file
rule merge_all_dnm_files:
    input:
        raw_denovos = "downsampling/csv/prefiltered/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.denovo.csv",
        phasing = "downsampling/csv/phasing/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.allele_sequences.tsv",
        transmission = "downsampling/csv/transmission/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.tsv",
        grandparental = "downsampling/csv/grandparents/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.tsv"
    output: out = "downsampling/csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.{KID_TARGET}.{MOM_TARGET}.{DAD_TARGET}.tsv"
    script:
        "rules/scripts/merge_mutations_with_metadata.py"
