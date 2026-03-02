import pandas as pd

# create global dictionaries we'll use
PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str})

CHILDREN = ped[(ped["paternal_id"] != "missing") & (~ped["paternal_id"].isin(["2281", "2214"]))]["sample_id"].to_list()
ALL_SAMPLES = ped["sample_id"].to_list()

ASSEMBLIES = ["CHM13v2"]
USE_NEW_BAMS = ["TOPUP", "ORIGINAL"]

wildcard_constraints:
    SAMPLE = r"[0-9]{4,6}",
    ASSEMBLY = "GRCh38|CHM13v2",
    TECH = "element|hifi|illumina",
    TOPUP = "TOPUP|ORIGINAL"

rule all:
    input:
        expand("csv/filtered_for_plots/{ASSEMBLY}.{USE_NEW_BAM}.tsv", ASSEMBLY=ASSEMBLIES, USE_NEW_BAM=USE_NEW_BAMS),
        expand("plots/{plot_name}.{ASSEMBLY}.{USE_NEW_BAM}.png", plot_name=[
            "length_and_purity", 
            "het_effect", 
            "rate_vs_age", 
            "parsimony_vs_inferred", 
            "dnm_counts", 
            "censat_rate", 
            "motif_sizes",
            "motif_rate",
            ], ASSEMBLY=ASSEMBLIES, USE_NEW_BAM=USE_NEW_BAMS),
        expand("plots/vntrs/{ASSEMBLY}.{TECH}.{USE_NEW_BAM}.plotted_trids.txt", ASSEMBLY=ASSEMBLIES, TECH=["hifi"], USE_NEW_BAM=USE_NEW_BAMS),
        expand("plots/svs/{ASSEMBLY}.{TECH}.{USE_NEW_BAM}.plotted_trids.txt", ASSEMBLY=ASSEMBLIES, TECH=["hifi"], USE_NEW_BAM=USE_NEW_BAMS),
        expand("plots/censat/{ASSEMBLY}.{TECH}.{USE_NEW_BAM}.plotted_trids.txt", ASSEMBLY=ASSEMBLIES, TECH=["hifi"], USE_NEW_BAM=USE_NEW_BAMS),
        expand("plots/recurrent/{ASSEMBLY}.{TECH}.{USE_NEW_BAM}.plotted_trids.txt", ASSEMBLY=ASSEMBLIES, TECH=["hifi"], USE_NEW_BAM=USE_NEW_BAMS),
        expand("plots/orthogonal_validation.{ASSEMBLY}.{TECH}.{USE_NEW_BAM}.png", ASSEMBLY=ASSEMBLIES, TECH=["element"], USE_NEW_BAM=USE_NEW_BAMS)



rule filter_denovos:
    input:
        dnms = expand("csv/annotated/{SAMPLE}.{{ASSEMBLY}}.element.{{TOPUP}}.tsv", SAMPLE=CHILDREN),
        denominators = expand("csv/denominators/{SAMPLE}.{{ASSEMBLY}}.{{TOPUP}}.denominator.tsv", SAMPLE=CHILDREN),
        recurrent = "csv/recurrent/{ASSEMBLY}.{TOPUP}.recurrent.tsv",
        censat = "data/t2t.censat.bed",
        fail = "analysis_and_plotting/FAILING_TRIDS.py"
    output:
        csv = "csv/filtered_for_plots/{ASSEMBLY}.{TOPUP}.tsv"
    params:
        filter_by_transmission = False,
        filter_by_grandparents = True,
        filter_by_orthogonal = True,
        filter_by_recurrent = True,
    script:
        "analysis_and_plotting/filter_dnms.py"



rule plot_dnm_counts:
    input:
        mutations = "csv/filtered_for_plots/{ASSEMBLY}.{TOPUP}.tsv"
    output:
        png = "plots/dnm_counts.{ASSEMBLY}.{TOPUP}.png"
    script:
        "analysis_and_plotting/plot_dnms.py"


rule plot_length_and_purity:
    input:
        mutations = "csv/filtered_for_plots/{ASSEMBLY}.{TOPUP}.tsv"
    output:
        png = "plots/length_and_purity.{ASSEMBLY}.{TOPUP}.png"
    script:
        "analysis_and_plotting/plot_purity.py"


rule plot_rate_vs_age:
    input:
        mutations = "csv/filtered_for_plots/{ASSEMBLY}.{TOPUP}.tsv",
        metadata = "tr_validation/data/k20_parental_age_at_birth.csv"
    output:
        png = "plots/rate_vs_age.{ASSEMBLY}.{TOPUP}.png"
    script:
        "analysis_and_plotting/plot_rate_vs_age.py"


rule plot_parsimony_vs_inferred:
    input:
        mutations = "csv/filtered_for_plots/{ASSEMBLY}.{TOPUP}.tsv",
    output:
        png = "plots/parsimony_vs_inferred.{ASSEMBLY}.{TOPUP}.png"
    script:
        "analysis_and_plotting/plot_parsimony_vs_inferred.py"


rule plot_het_effect:
    input:
        mutations = "csv/filtered_for_plots/{ASSEMBLY}.{TOPUP}.tsv",
    output:
        png = "plots/het_effect.{ASSEMBLY}.{TOPUP}.png"
    script:
        "analysis_and_plotting/plot_het_effect.py"


rule plot_motif_sizes:
    input:
        mutations = "csv/filtered_for_plots/{ASSEMBLY}.{TOPUP}.tsv",
        akshay = "data/K1463.CHM13v2.DNMs.416.demintr.output",
    output:
        by_motif = "plots/motif_sizes.{ASSEMBLY}.{TOPUP}.png",
        by_tr_type = "plots/tr_type_sizes.{ASSEMBLY}.{TOPUP}.png"
    script:
        "analysis_and_plotting/plot_dnm_motif_sizes.py"


rule plot_censat_rate:
    input:
        mutations = "csv/filtered_for_plots/{ASSEMBLY}.{TOPUP}.tsv",
        denominators = expand("csv/denominators/{SAMPLE}.{{ASSEMBLY}}.{{TOPUP}}.denominator.tsv", SAMPLE=CHILDREN),
    output:
        png = "plots/censat_rate.{ASSEMBLY}.{TOPUP}.png"
    params:
        rate_per_haplotype = True
    script:
        "analysis_and_plotting/plot_rate_in_censat.py"


rule plot_rate_by_motif:
    input:
        mutations = "csv/filtered_for_plots/{ASSEMBLY}.{TOPUP}.tsv",
        denominators = expand("csv/denominators/{SAMPLE}.{{ASSEMBLY}}.{{TOPUP}}.denominator.tsv", SAMPLE=CHILDREN),
    output:
        png = "plots/motif_rate.{ASSEMBLY}.{TOPUP}.png"
    params:
        rate_per_haplotype = True
    script:
        "analysis_and_plotting/plot_rate_vs_motif.py"


rule plot_read_evidence_vntr:
    input:
        mutations = expand("csv/annotated/{SAMPLE}.{{ASSEMBLY}}.{{TECH}}.{{TOPUP}}.tsv", SAMPLE=CHILDREN),
        recurrent = "csv/recurrent/{ASSEMBLY}.{TOPUP}.recurrent.tsv",
                censat = "data/t2t.censat.bed"

    output:
        fh = "plots/vntrs/{ASSEMBLY}.{TECH}.{TOPUP}.plotted_trids.txt"
    params:
        minimum_motif_size = 7,
        minimum_dnm_size = 1,
        outpref = "plots/vntrs",
                censat_only = False,

    script:
        "analysis_and_plotting/plot_read_evidence.py"


rule plot_read_evidence_svs:
    input:
        mutations = expand("csv/annotated/{SAMPLE}.{{ASSEMBLY}}.{{TECH}}.{{TOPUP}}.tsv", SAMPLE=CHILDREN),
        recurrent = "csv/recurrent/{ASSEMBLY}.{TOPUP}.recurrent.tsv",
                censat = "data/t2t.censat.bed"

    output:
        fh = "plots/svs/{ASSEMBLY}.{TECH}.{TOPUP}.plotted_trids.txt"
    params:
        minimum_motif_size = 1,
        minimum_dnm_size = 50,
        outpref = "plots/svs",
                censat_only = False,

    script:
        "analysis_and_plotting/plot_read_evidence.py"


rule plot_read_evidence_censat:
    input:
        mutations = expand("csv/annotated/{SAMPLE}.{{ASSEMBLY}}.{{TECH}}.{{TOPUP}}.tsv", SAMPLE=CHILDREN),
        recurrent = "csv/recurrent/{ASSEMBLY}.{TOPUP}.recurrent.tsv",
        censat = "data/t2t.censat.bed"
    output:
        fh = "plots/censat/{ASSEMBLY}.{TECH}.{TOPUP}.plotted_trids.txt"
    params:
        minimum_motif_size = 1,
        minimum_dnm_size = 1,
        outpref = "plots/censat",
        censat_only = True,
    script:
        "analysis_and_plotting/plot_read_evidence.py"


rule plot_read_evidence_recurrent:
    input:
        mutations = expand("csv/orthogonal_support/all/{SAMPLE}.{{ASSEMBLY}}.{{TECH}}.{{TOPUP}}.tsv", SAMPLE=ALL_SAMPLES),
        recurrent = "csv/recurrent/{ASSEMBLY}.{TOPUP}.recurrent.tsv",
        ped = PED_FILE
    output:
        fh = "plots/recurrent/{ASSEMBLY}.{TECH}.{TOPUP}.plotted_trids.txt"
    params:
        outpref = "plots/recurrent"
    script:
        "analysis_and_plotting/plot_read_evidence_recurrent.py"


rule calculate_orthogonal_validation:
    input:
        mutations = expand("csv/annotated/{SAMPLE}.{{ASSEMBLY}}.{{TECH}}.{{TOPUP}}.tsv", SAMPLE=CHILDREN),
        recurrent = "csv/recurrent/{ASSEMBLY}.{TOPUP}.recurrent.tsv",
        censat = "data/t2t.censat.bed"
    output:
        png = "plots/orthogonal_validation.{ASSEMBLY}.{TECH}.{TOPUP}.png"
    params:
        minimum_motif_size = 7
    script:
        "analysis_and_plotting/calculate_orthogonal_validation.py"
        
