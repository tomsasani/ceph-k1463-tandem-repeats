import pandas as pd

# create global dictionaries we'll use
PED_FILE = "tr_validation/data/file_mapping.csv"
ped = pd.read_csv(PED_FILE, sep=",", dtype={"paternal_id": str, "maternal_id": str, "sample_id": str})

CHILDREN = ped[(ped["paternal_id"] != "missing") & (~ped["paternal_id"].isin(["2281", "2214"]))]["sample_id"].to_list()

ASSEMBLIES = ["CHM13v2", "GRCh38"]

wildcard_constraints:
    SAMPLE = r"[0-9]{4,6}",
    ASSEMBLY = "GRCh38|CHM13v2",
    TECH = "element|hifi|illumina",
    TOPUP = "TOPUP|ORIGINAL"

rule all:
    input:
        expand("csv/filtered_for_plots/{ASSEMBLY}.TOPUP.tsv", ASSEMBLY=ASSEMBLIES),
        expand("plots/{plot_name}.{ASSEMBLY}.TOPUP.png", plot_name=[
            "length_and_purity", 
            "het_effect", 
            "rate_vs_age", 
            "parsimony_vs_inferred", 
            "dnm_counts", 
            "censat_rate", 
            "motif_sizes",
            "motif_rate",
            ], ASSEMBLY=ASSEMBLIES)


rule filter_denovos:
    input:
        dnms = expand("csv/annotated/{SAMPLE}.{{ASSEMBLY}}.element.{{TOPUP}}.tsv", SAMPLE=CHILDREN),
        denominators = expand("csv/denominators/{SAMPLE}.{{ASSEMBLY}}.{{TOPUP}}.denominator.tsv", SAMPLE=CHILDREN),
        recurrent = "csv/recurrent/{ASSEMBLY}.{TOPUP}.recurrent.tsv",
        censat = "data/t2t.censat.bed"
    output:
        csv = "csv/filtered_for_plots/{ASSEMBLY}.{TOPUP}.tsv"
    params:
        filter_by_transmission = False,
        filter_by_grandparents = False,
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
        png = "plots/motif_sizes.{ASSEMBLY}.{TOPUP}.png"
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
