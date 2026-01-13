CHROMS = list(map(str, range(1, 23)))
CHROMS = [f"chr{c}" for c in CHROMS]
CHROMS.extend(["chrX", "chrY"])

# create global dictionaries we'll use
PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str})


ALL_SAMPLES = ped["sample_id"].to_list()

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


rule call_snvs:
    input:
        ref = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz",
        ref_idx = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz.fai",
        bam = lambda wildcards: get_complete_bams(wildcards),
        bam_idx = lambda wildcards: get_complete_bams(wildcards) + ".bai",
        sif = "deepvariant_1.9.0.sif"
    output:
        vcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.vcf.gz",
        gvcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.g.vcf.gz",
    threads: 8
    params:
        alt_sample_id = lambda wildcards: SMP2ALT[wildcards.SAMPLE]
    resources:
        mem_mb = 64_000,
    script:
        "bash_scripts/run_deepvariant.sh"


rule sort_snv_chrom_vcfs:
    input: 
        vcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.vcf.gz",
    output: "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz"
    threads: 4
    shell:
        """
        module load bcftools

        bcftools sort -Ob -o {output} {input.vcf}
        """


rule index_snv_chrom_vcfs:
    input: vcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz"
    output: "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz.tbi"
    threads: 4
    script:
        "bash_scripts/index_vcf.sh"


rule combine_trio_chrom_vcfs:
    input:
        sif = "glnexus_v1.2.7.sif",
        kid_snp_gvcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.g.vcf.gz",
        dad_snp_gvcf = lambda wildcards: f"data/vcf/snv/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.g.vcf.gz",
        mom_snp_gvcf = lambda wildcards: f"data/vcf/snv/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.g.vcf.gz",
    output: "data/vcf/trios/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.trio.vcf.gz"
    threads: 4
    params:
        gl_nexus_prefix = lambda wildcards: f"gl_nexus_dbs/{wildcards.SAMPLE}_{wildcards.ASSEMBLY}_{wildcards.USE_NEW_BAM}_{wildcards.CHROM}"
    resources:
        mem_mb = 32_000
    script:
        "bash_scripts/combine_trio_vcf.sh"


rule combine_cohort_chrom_vcfs:
    input:
        sif = "glnexus_v1.2.7.sif",
        all_gvcfs = expand("data/vcf/snv/per-chrom/{SAMPLE}.{{ASSEMBLY}}.{{USE_NEW_BAM}}.{{CHROM}}.g.vcf.gz", SAMPLE=ALL_SAMPLES),
    output: "data/vcf/cohort/{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.cohort.vcf.gz"
    threads: 4
    params:
        gl_nexus_prefix = lambda wildcards: f"gl_nexus_dbs/cohort_{wildcards.ASSEMBLY}_{wildcards.USE_NEW_BAM}_{wildcards.CHROM}"
    resources:
        mem_mb = 32_000
    script:
        "bash_scripts/combine_cohort_vcf.sh"


rule combine_cohort_vcfs:
    input:
        vcfs = expand("data/vcf/cohort/{{ASSEMBLY}}.{{USE_NEW_BAM}}.{CHROM}.cohort.vcf.gz", CHROM=CHROMS),
    output: "data/vcf/snv/{ASSEMBLY}.{USE_NEW_BAM}.cohort.vcf.gz"
    resources:
        mem_mb = 32_000
    threads: 8
    shell:
        """
        module load bcftools
        
        bcftools concat {input.vcfs} | bcftools view --threads {threads} | bgzip > {output}
        """


rule merge_trio_vcfs:
    input:
        vcfs = expand("data/vcf/trios/per-chrom/{{SAMPLE}}.{{ASSEMBLY}}.{{USE_NEW_BAM}}.{CHROM}.trio.vcf.gz", CHROM=CHROMS)
    output: "data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.vcf.gz"
    threads: 8
    shell:
        """
        module load bcftools
        
        bcftools concat {input.vcfs} | bcftools view --threads {threads} | bgzip > {output}
        """


rule index_merged_vcf:
    input: vcf = "data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.vcf.gz"
    output: "data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.vcf.gz.tbi"
    script:
        "bash_scripts/index_vcf.sh"