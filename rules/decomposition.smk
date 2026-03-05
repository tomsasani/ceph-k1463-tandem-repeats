def get_hprc_sample_vcf(wildcards):
    assembly_adj = wildcards.ASSEMBLY.lower() if "CHM" in wildcards.ASSEMBLY else wildcards.ASSEMBLY
    return f"{HPRC_PREF}/{wildcards.SAMPLE}_{assembly_adj}_50bp_merge.sorted.vcf.gz"

def get_ceph_sample_vcf(wildcards):
    return f"data/trgt/phased/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.TOPUP.phased.vcf.gz"

def get_ceph_sample_dnms(wildcards):
    if wildcards.SAMPLE in CHILDREN:
        return f"csv/prefiltered/denovos/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.TOPUP.tsv" 
    else:
        return None

rule extract_ceph_fasta:
    input:
        recurrent_dnms = "csv/recurrent/{ASSEMBLY}.TOPUP.recurrent.tsv",
        sample_vcf = get_ceph_sample_vcf,
        sample_dnms = get_ceph_sample_dnms
    output:
        fasta = temp("data/recurrent_fasta/ceph/{TRID}.{SAMPLE}.{ASSEMBLY}.fa")
    params:
    script:
        "scripts/extract_fasta_for_trviz.py"
        

rule extract_hprc_fasta:
    input:
        recurrent_dnms = "csv/recurrent/{ASSEMBLY}.TOPUP.recurrent.tsv",
        sample_vcf = get_hprc_sample_vcf,
    output:
        fasta = temp("data/recurrent_fasta/hprc/{TRID}.{SAMPLE}.{ASSEMBLY}.fa")
    params:
        use_dnms = False
    script:
        "scripts/extract_fasta_for_trviz.py"


rule combine_ceph_fasta:
    input:
        fastas = expand("data/recurrent_fasta/ceph/{{TRID}}.{SAMPLE}.{{ASSEMBLY}}.fa", SAMPLE=ALL_SAMPLES)
    output:
        fasta = "data/recurrent_fasta/ceph/combined/{TRID}.{ASSEMBLY}.combined.fa"
    shell:
        """
        cat {input.fastas} > {output.fasta}
        """

rule combine_hprc_fasta:
    input:
        fastas = expand("data/recurrent_fasta/hprc/{{TRID}}.{SAMPLE}.{{ASSEMBLY}}.fa", SAMPLE=HPRC_SAMPLES)
    output:
        fasta = "data/recurrent_fasta/hprc/combined/{TRID}.{ASSEMBLY}.combined.fa"
    shell:
        """
        cat {input.fastas} > {output.fasta}
        """

rule decompose:
    input:
        recurrent_dnms = "csv/recurrent/{ASSEMBLY}.TOPUP.recurrent.tsv",
        fasta = "data/recurrent_fasta/{COHORT}/combined/{TRID}.{ASSEMBLY}.combined.fa"
    output:
        png = "trviz/{COHORT}/{TRID}.{ASSEMBLY}.png",
        seq_tsv = "trviz/{COHORT}/{TRID}.{ASSEMBLY}.encoded.tsv",
        key_tsv = "trviz/{COHORT}/{TRID}.{ASSEMBLY}.key.tsv",
        aln_fa = "trviz/{COHORT}/{TRID}.{ASSEMBLY}_alignment_output.fa"
    script:
        "scripts/decompose_trviz.py"