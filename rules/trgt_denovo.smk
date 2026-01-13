CUR_PREF = "/scratch/ucgd/lustre-labs/quinlan/u1006375/CEPH-K1463-TandemRepeats/"

rule run_trgt_denovo:
    input:
        reference = CUR_PREF + "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz",
        repeats =  CUR_PREF + "data/catalogs/{CHROM}.{ASSEMBLY}.catalog.bed",
        trgt_denovo_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/trgt-denovo",
        kid_vcf = CUR_PREF + "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz",
        kid_vcf_idx = CUR_PREF + "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz.tbi",
        mom_vcf = lambda wildcards: CUR_PREF + f"data/trgt/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.sorted.vcf.gz",
        mom_vcf_idx = lambda wildcards: CUR_PREF + f"data/trgt/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.sorted.vcf.gz.tbi",
        dad_vcf = lambda wildcards: CUR_PREF + f"data/trgt/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.sorted.vcf.gz",
        dad_vcf_idx = lambda wildcards: CUR_PREF + f"data/trgt/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.sorted.vcf.gz.tbi",
        kid_bam = CUR_PREF + "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.spanning.sorted.bam",
        kid_bam_idx = CUR_PREF + "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.spanning.sorted.bam.bai",
        mom_bam = lambda wildcards: CUR_PREF + f"data/trgt/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.spanning.sorted.bam",
        mom_bam_idx = lambda wildcards: CUR_PREF + f"data/trgt/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.spanning.sorted.bam.bai",
        dad_bam = lambda wildcards: CUR_PREF + f"data/trgt/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.spanning.sorted.bam",
        dad_bam_idx = lambda wildcards: CUR_PREF + f"data/trgt/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.spanning.sorted.bam.bai",
    output:
        output = CUR_PREF + "csv/raw_denovos/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.trgt-denovo.csv"
    threads: 32
    resources:
        mem_mb = 32_000,
        runtime = 720,
    params:
        kid_pref = lambda wildcards: CUR_PREF + f"data/trgt/per-chrom/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}",
        mom_pref = lambda wildcards: CUR_PREF + f"data/trgt/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}",
        dad_pref = lambda wildcards: CUR_PREF + f"data/trgt/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}",
        output_dir = lambda wildcards: f"trgt_denovo_tmp/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}"
    script:
        "bash_scripts/run_trgt_denovo.sh"