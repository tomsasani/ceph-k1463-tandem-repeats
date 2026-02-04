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
        vcf = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz",
    script:
        "bash_scripts/sort_trgt_vcf.sh"


rule index_trgt_vcfs:
    input:
        vcf = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz"
    output:
        "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz.tbi"
    script:
        "bash_scripts/index_vcf.sh"


rule sort_bam:
    input: bam = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.spanning.bam"
    output: bam = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.spanning.sorted.bam"
    script:
        "bash_scripts/sort_bam.sh"


rule index_bam:
    input: bam = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.spanning.sorted.bam"
    output: "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.spanning.sorted.bam.bai"
    script:
        "bash_scripts/index_bam.sh"