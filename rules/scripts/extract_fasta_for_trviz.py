from cyvcf2 import VCF
import pandas as pd
import numpy as np

recurrent = pd.read_csv(snakemake.input.recurrent_dnms, sep="\t")
recurrent = recurrent[recurrent["trid"] == snakemake.wildcards.TRID]

# figure out which haplotype has the DNM if this is a kid
denovo_gt = None
if snakemake.params.sample_dnms is not None:
    sample_dnms = pd.read_csv(snakemake.params.sample_dnms, sep="\t")
    dnms = sample_dnms[sample_dnms["trid"] == snakemake.wildcards.TRID]
    # get index with de novo
    denovo_gt = dnms["genotype"].values

vcf = VCF(snakemake.params.sample_vcf, gts012=True)

outfh = open(snakemake.output.fasta, "w")

for i, row in recurrent.iterrows():
    chrom, start, end = row["#chrom"], row["start"], row["end"]
    region = f"{chrom}:{start}-{end}"
    for v in vcf(region):
        if v.gt_types[0] == 3: continue
        hap_a, hap_b = np.array(v.genotypes)[0, :-1]
        ref, alts = v.REF, v.ALT
        alleles = [ref] + alts

        for hap_i, hap_gt in enumerate((hap_a, hap_b)):
            denovo_status = "nondenovo"
            if denovo_gt is not None and hap_gt in denovo_gt: 
                denovo_status = "denovo"
            header = f">{snakemake.wildcards.SAMPLE}_{snakemake.wildcards.TRID}_haplotype_{hap_i}_{denovo_status}"
            print ("\n".join([header, alleles[hap_gt]]), file=outfh)

        
outfh.close()