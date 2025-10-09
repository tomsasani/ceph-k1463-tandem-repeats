rule run_trgt_denovo:
    input:
        reference = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz",
        repeats =  "data/catalogs/{CHROM}.{ASSEMBLY}.catalog.bed",
        trgt_denovo_binary = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/trgt-denovo",
        kid_vcf = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz",
        kid_vcf_idx = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz.tbi",
        mom_vcf = lambda wildcards: f"data/trgt/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.sorted.vcf.gz",
        mom_vcf_idx = lambda wildcards: f"data/trgt/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.sorted.vcf.gz.tbi",
        dad_vcf = lambda wildcards: f"data/trgt/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.sorted.vcf.gz",
        dad_vcf_idx = lambda wildcards: f"data/trgt/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.sorted.vcf.gz.tbi",
        kid_bam = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.spanning.sorted.bam",
        mom_bam = lambda wildcards: f"data/trgt/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.spanning.sorted.bam",
        dad_bam = lambda wildcards: f"data/trgt/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}.spanning.sorted.bam",
    output:
        output = "csv/raw_denovos/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.trgt-denovo.csv"
    threads: 32
    resources:
        mem_mb = 32_000,
        runtime = 720,
        slurm_account = "ucgd-rw",
        slurm_partition = "ucgd-rw"
    params:
        kid_pref = lambda wildcards: f"data/trgt/per-chrom/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}",
        mom_pref = lambda wildcards: f"data/trgt/per-chrom/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}",
        dad_pref = lambda wildcards: f"data/trgt/per-chrom/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.{wildcards.CHROM}",
    script:
        "bash_scripts/run_trgt_denovo.sh"