rule run_hiphase:
    input:
        snv_vcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz",
        snv_vcf_idx = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz.tbi",
        str_vcf = CUR_PREF + "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz",
        str_vcf_idx = CUR_PREF + "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.sorted.vcf.gz.tbi",
        bam = lambda wildcards: get_complete_bams(wildcards),
        bam_idx = lambda wildcards: get_complete_bams(wildcards) + ".bai",
        reference = "data/contigs/{CHROM}.{ASSEMBLY}.fa.gz",
    output:
        snv_vcf = "data/vcf/snv/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.phased.vcf.gz",
        str_vcf = "data/trgt/per-chrom/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.{CHROM}.phased.vcf.gz",
    threads: 4
    resources:
        mem_mb = 32_000
    shell:
        """
        ~/bin/hiphase-v1.4.4-x86_64-unknown-linux-gnu/hiphase --threads {threads} \
                      --bam {input.bam} \
                      --vcf {input.snv_vcf} \
                      --vcf {input.str_vcf} \
                      --reference {input.reference} \
                      --output-vcf {output.snv_vcf} \
                      --output-vcf {output.str_vcf} \

        """


rule combine_phased_snv_vcfs:
    input:
        vcfs = expand("data/vcf/snv/per-chrom/{{SAMPLE}}.{{ASSEMBLY}}.{{USE_NEW_BAM}}.{CHROM}.phased.vcf.gz", CHROM=CHROMS)
    output: "data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.vcf.gz"
    threads: 8
    shell:
        """
        module load bcftools
        
        bcftools concat {input.vcfs} | bcftools view --threads {threads} -f 'PASS' | bgzip > {output}
        """


rule index_phased_snv_vcfs:
    input: vcf = "data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.vcf.gz"
    output: "data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.vcf.gz.tbi"
    script:
        "bash_scripts/index_vcf.sh"


rule combine_phased_trgt_vcfs:
    input:
        vcfs = expand("data/trgt/per-chrom/{{SAMPLE}}.{{ASSEMBLY}}.{{USE_NEW_BAM}}.{CHROM}.phased.vcf.gz", CHROM=CHROMS)
    output: "data/trgt/phased/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.vcf.gz"
    threads: 8
    shell:
        """
        module load bcftools
        
        bcftools concat {input.vcfs} | bcftools view --threads {threads} | bgzip > {output}
        """


rule index_phased_trgt_vcfs:
    input: vcf = "data/trgt/phased/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.vcf.gz"
    output: "data/trgt/phased/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.vcf.gz.tbi"
    script:
       "bash_scripts/index_vcf.sh"