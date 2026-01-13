import pandas as pd
import glob

CUR_PREF = "/scratch/ucgd/lustre-labs/quinlan/u1006375/CEPH-K1463-TandemRepeats/"

wildcard_constraints:
    SAMPLE = r"[0-9]{4,6}",
    ASSEMBLY = "GRCh38|CHM13v2",
    TECH = r"[a-z]+",
    USE_NEW_BAM = "TOPUP|ORIGINAL"


# create global dictionaries we'll use
PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str})


SMP2DAD = dict(zip(ped["sample_id"], ped["paternal_id"]))
SMP2MOM = dict(zip(ped["sample_id"], ped["maternal_id"]))
SMP2ALT = dict(zip(ped["sample_id"], ped["alt_sample_id"]))
SMP2SUFF = dict(zip(ped["sample_id"], ped["suffix"]))


CHILDREN = ped[(ped["paternal_id"] != "missing") & (~ped["paternal_id"].isin(["2281", "2214"]))]["sample_id"].to_list()
ALL_SAMPLES = ped["sample_id"].to_list()


HPRC_PREF = f"/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/"

HPRC_SAMPLES = []
for fh in glob.glob(HPRC_PREF + "GRCh38_v1.0_50bp_merge/1.1.2-69937d83/hprc/*.vcf.gz"):
    sample_id = fh.split("/")[-1].split("_")[0]
    if sample_id == "hprc": continue
    HPRC_SAMPLES.append(sample_id)


def get_bam_for_validation(sample, assembly, tech):

    assembly_adj = assembly.split('v2')[0]
    sample_id = SMP2ALT[sample] if tech in ("ont", "hifi") else sample

    tech2path = {
        "ont": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/ont-bams/{0}/{1}.minimap2.bam".format(assembly_adj, sample_id,), 
        "element": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/element/{0}_bams/{1}-E.{0}.merged.sort.bam".format(assembly_adj, sample_id,),
        "illumina": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/CEPH/cram/{0}.cram".format(sample_id),
        "hifi": "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams/{0}/{1}.{2}.haplotagged.bam".format(assembly_adj, sample_id, assembly_adj.lower() if 'CHM' in assembly else assembly_adj,),
        }

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
    if wildcards.ASSEMBLY == "GRCh38":
        return "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/human_GRCh38_no_alt_analysis_set.palladium-v1.0.trgt.annotations.bed.gz"
        #return "data/gangstr.annotations.bed.gz"
    elif wildcards.ASSEMBLY == "CHM13v2":
        return "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/TRGT/from_aws/staging/catalogs/chm13v2.0_maskedY_rCRS.palladium-v1.0.trgt.annotations.bed.gz"


rule combine_trgt_denovo:
    input:
        fhs = expand(CUR_PREF + "csv/raw_denovos/{{SAMPLE}}.{{ASSEMBLY}}.{{USE_NEW_BAM}}.{CHROM}.trgt-denovo.csv", CHROM=CHROMS)
    output: fh = "csv/denovos/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.trgt-denovo.csv"
    shell:
        """
        grep '^chrom' {input.fhs[0]} > {output.fh}

        cat {input.fhs} | grep -v '^chrom'  >> {output.fh}
        """


rule prefilter_all_loci:
    input:
        kid_mutation_df = "csv/denovos/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.trgt-denovo.csv",
        ped = "tr_validation/data/file_mapping.csv",
        annotations = lambda wildcards: get_annotation_fh(wildcards),
        utils = "rules/scripts/utils.py"

    output: 
        fh = "csv/prefiltered/all_loci/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.tsv"
    params:
        filtering_mode = "lenient"
    resources:
        mem_mb = 32_000
    script:
        "scripts/filter_mutation_dataframe.py"


rule prefilter_denovos:
    input:
        ped = "tr_validation/data/file_mapping.csv",
        kid_mutation_df = "csv/denovos/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.trgt-denovo.csv",
        annotations = lambda wildcards: get_annotation_fh(wildcards),
        utils = "rules/scripts/utils.py"
    output: 
        fh = "csv/prefiltered/denovos/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.tsv"
    params:
        filtering_mode = "strict"
    resources:
        mem_mb = 32_000
    script:
        "scripts/filter_mutation_dataframe.py"


rule calculate_denominator:
    input:
        loci = "csv/prefiltered/all_loci/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.tsv",
        censat = "data/t2t.censat.bed"
    output:  out = "csv/denominators/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.denominator.tsv"
    script:
        "scripts/calculate_denominator.py"


rule annotate_with_informative_sites:
    input:
        cohort_snp_vcf = "data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.vcf.gz",
        cohort_snp_vcf_idx = "data/vcf/trios/merged/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.vcf.gz.tbi",
        kid_phased_snp_vcf = "data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.vcf.gz",
        kid_phased_snp_vcf_idx = "data/vcf/snv/phased/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.vcf.gz.tbi",
        kid_mutation_df = "csv/prefiltered/denovos/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.tsv",
        kid_phased_str_vcf = "data/trgt/phased/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.vcf.gz",
        kid_phased_str_vcf_idx = "data/trgt/phased/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.vcf.gz.tbi",
        py_script = "rules/scripts/annotate_with_informative_sites.py",
    output:
        out = "csv/phasing/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.annotated.2gen.tsv"
    params:
        dad_id = lambda wildcards: SMP2ALT[SMP2DAD[wildcards.SAMPLE]],
        mom_id = lambda wildcards: SMP2ALT[SMP2MOM[wildcards.SAMPLE]],
        focal_alt_id = lambda wildcards: SMP2ALT[wildcards.SAMPLE],
        slop = 500_000
    script:
        "scripts/annotate_with_informative_sites.py"
        


rule annotate_with_parental_haplotype:
    input:
        annotated_dnms = "csv/phasing/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.annotated.2gen.tsv",
        dad_phased_str_vcf = lambda wildcards: f"data/trgt/phased/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.phased.vcf.gz",
        dad_phased_str_vcf_idx = lambda wildcards: f"data/trgt/phased/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.phased.vcf.gz.tbi",
        dad_phased_snv_vcf = lambda wildcards: f"data/vcf/snv/phased/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.phased.vcf.gz",
        dad_phased_snv_vcf_idx = lambda wildcards: f"data/vcf/snv/phased/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.phased.vcf.gz.tbi",
        mom_phased_str_vcf = lambda wildcards: f"data/trgt/phased/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.phased.vcf.gz",
        mom_phased_str_vcf_idx = lambda wildcards: f"data/trgt/phased/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.phased.vcf.gz.tbi",
        mom_phased_snv_vcf = lambda wildcards: f"data/vcf/snv/phased/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.phased.vcf.gz",
        mom_phased_snv_vcf_idx = lambda wildcards: f"data/vcf/snv/phased/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.phased.vcf.gz.tbi",
    output:
        out = "csv/phasing/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.haplotyped.2gen.tsv"
    params:
        dad_id = lambda wildcards: SMP2ALT[SMP2DAD[wildcards.SAMPLE]],
        mom_id = lambda wildcards: SMP2ALT[SMP2MOM[wildcards.SAMPLE]],
    script:
        "scripts/annotate_with_parental_haplotype.py"


rule phase:
    input:
        annotated_dnms = "csv/phasing/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.haplotyped.2gen.tsv",
        py_script = "rules/scripts/phase_by_ps.py",
    output: "csv/phasing/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.2gen.tsv"
    shell:
        """
        python {input.py_script} --annotated_dnms {input.annotated_dnms} \
                                 --out {output}
        """


rule phase_three_gen:
    input:
        mutation_df = "csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.tsv",
        cohort_vcf = "data/vcf/snv/CHM13v2.TOPUP.cohort.vcf.gz"
    output:
        tsv = "csv/phasing/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.3gen.tsv"
    params:
        slop = 250_000
    script:
        "scripts/phase_three_gen.py"


rule annotate_with_parental_alleles:
    input:
        annotated_dnms = "csv/phasing/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.phased.2gen.tsv",
        dad_phased_str_vcf = lambda wildcards: f"data/trgt/phased/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.phased.vcf.gz",
        dad_phased_str_vcf_idx = lambda wildcards: f"data/trgt/phased/{SMP2DAD[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.phased.vcf.gz.tbi",
        mom_phased_str_vcf = lambda wildcards: f"data/trgt/phased/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.phased.vcf.gz",
        mom_phased_str_vc_idx = lambda wildcards: f"data/trgt/phased/{SMP2MOM[wildcards.SAMPLE]}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.phased.vcf.gz.tbi",
        kid_phased_str_vcf = lambda wildcards: f"data/trgt/phased/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.phased.vcf.gz",
        kid_phased_str_vcf_idx = lambda wildcards: f"data/trgt/phased/{wildcards.SAMPLE}.{wildcards.ASSEMBLY}.{wildcards.USE_NEW_BAM}.phased.vcf.gz.tbi",
    output:
        out = "csv/phasing/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.allele_sequences.tsv"
    script:
        "scripts/annotate_with_parental_allele_sequences.py"



rule validate_dnms_with_orthogonal_tech:
    input:
        mutations = "csv/prefiltered/denovos/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.tsv"
    output: fh = "csv/orthogonal_support/{SAMPLE}.{ASSEMBLY}.{TECH}.{USE_NEW_BAM}.tsv",
    params:
        kid_bam = lambda wildcards: get_bam_for_validation(wildcards.SAMPLE, wildcards.ASSEMBLY, wildcards.TECH),
        mom_bam = lambda wildcards: get_bam_for_validation(SMP2MOM[wildcards.SAMPLE], wildcards.ASSEMBLY, wildcards.TECH),
        dad_bam = lambda wildcards: get_bam_for_validation(SMP2DAD[wildcards.SAMPLE], wildcards.ASSEMBLY, wildcards.TECH),
    script:
        "scripts/annotate_with_orthogonal_evidence.py"


rule combine_raw_dnms:
    input:
        kid_mutation_dfs = expand("csv/prefiltered/denovos/{SAMPLE}.{{ASSEMBLY}}.{{USE_NEW_BAM}}.tsv", SAMPLE=CHILDREN)
    output:
        fh = "csv/all_dnms/{ASSEMBLY}.{USE_NEW_BAM}.tsv"
    run:
        import pandas as pd
        res = []
        for fh in input.kid_mutation_dfs:
            df = pd.read_csv(fh, sep="\t")
            res.append(df)
        res = pd.concat(res)
        res.to_csv(output.fh, sep="\t", index=False)
        

rule validate_all_with_orthogonal_tech:
    input:
        mutations = "csv/all_dnms/{ASSEMBLY}.{USE_NEW_BAM}.tsv"
    output: fh = "csv/orthogonal_support/all/{SAMPLE}.{ASSEMBLY}.{TECH}.{USE_NEW_BAM}.tsv",
    params:
        kid_bam = lambda wildcards: get_bam_for_validation(wildcards.SAMPLE, wildcards.ASSEMBLY, wildcards.TECH),
        mom_bam = None,
        dad_bam = None,
    script:
        "scripts/annotate_with_orthogonal_evidence.py"


rule add_orthogonal_filter:
    input:
        mutations = "csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.tsv",
        orthogonal_evidence = "csv/orthogonal_support/{SAMPLE}.{ASSEMBLY}.{TECH}.{USE_NEW_BAM}.tsv"
    output:
        out = "csv/annotated/{SAMPLE}.{ASSEMBLY}.{TECH}.{USE_NEW_BAM}.tsv"
    script:
        "scripts/add_orthogonal_filter.py"
    

rule add_transmission_evidence:
    input:
        mutations = "csv/prefiltered/denovos/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.tsv",
        other_vcfs = get_children_vcfs,
        other_vcf_idxs = lambda wildcards: [v + ".tbi" for v in get_children_vcfs(wildcards)]

    output: fh = "csv/transmission/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.tsv"
    params:
        generation_to_query = "children"
    script:
        "scripts/assess_presence_in_cohort.py"


rule add_grandparental_evidence:
    input:
        mutations = "csv/prefiltered/denovos/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.tsv",
        other_vcfs = lambda wildcards: get_grandparent_vcfs(wildcards),
        other_vcf_idxs = lambda wildcards: [v + ".tbi" for v in get_grandparent_vcfs(wildcards)]
    output: fh = "csv/grandparents/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.tsv"
    params:
        generation_to_query = "grandparents"
    script:
        "scripts/assess_presence_in_cohort.py"


rule merge_all_dnm_files:
    input:
        raw_denovos = "csv/prefiltered/denovos/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.tsv",
        phasing = "csv/phasing/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.allele_sequences.tsv",
        denominator = "csv/denominators/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.denominator.tsv",
        transmission = "csv/transmission/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.tsv",
        grandparental = "csv/grandparents/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.tsv",
        utils = "rules/scripts/utils.py"
    output: out = "csv/filtered_and_merged/{SAMPLE}.{ASSEMBLY}.{USE_NEW_BAM}.tsv"
    params: phase_threshold = 1.0
    script:
        "scripts/merge_mutations_with_metadata.py"


rule find_recurrents:
    input:
        fhs = expand("csv/filtered_and_merged/{SAMPLE}.{{ASSEMBLY}}.{{USE_NEW_BAM}}.tsv", SAMPLE=CHILDREN)
    output:
        out = "csv/recurrent/{ASSEMBLY}.{USE_NEW_BAM}.recurrent.tsv"
    script:
        "scripts/find_recurrent_trids.py"


rule annotate_with_hprc_heterozygosity:
    input:
        vcf = lambda wildcards: HPRC_PREF + f"{wildcards.ASSEMBLY}_v1.0_50bp_merge/1.1.2-69937d83/hprc/{wildcards.HPRC_SAMPLE}_{wildcards.ASSEMBLY}_50bp_merge.sorted.vcf.gz",
        denovos = "{ASSEMBLY}.filtered.tsv",
        catalog = lambda wildcards: ASSEMBLY2CATALOG[wildcards.ASSEMBLY],
    output:
        mutations = "data/hprc/{HPRC_SAMPLE}.{ASSEMBLY}.{KIND}.heterozygosity.tsv"     
    script:
        "scripts/calculate_hprc_heterozygosity.py"


rule merge_hprc:
    input:
        fhs = expand("data/hprc/{HPRC_SAMPLE}.{{ASSEMBLY}}.{KIND}.heterozygosity.tsv", HPRC_SAMPLE=HPRC_SAMPLES, KIND=["wt", "denovo"])
    output: fh = "csv/hprc/combined.{ASSEMBLY}.heterozygosity.tsv"
    shell:
        """
        cat {input.fhs} | grep -v 'chrom' > {output.fh}
        """