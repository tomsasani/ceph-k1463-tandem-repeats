import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import scipy.stats as ss
import tqdm

from snakemake.script import snakemake

plt.rc("font", size=12)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

mutations = pd.read_csv(snakemake.input.mutations, sep="\t", dtype={"paternal_id": str})

phased = mutations[mutations["phase"] != "unknown"]

hap_phased = phased[~phased["haplotype_in_parent_consensus"].str.contains("unknown")]

hap_phased["pass"] = hap_phased["haplotype_in_parent_consensus"].apply(lambda s: "unknown" not in s)
hap_phased = hap_phased[hap_phased["pass"] == True]

matching = hap_phased[hap_phased["likely_denovo_size"] == hap_phased["likely_denovo_size_parsimony"]]
matching["Both methods agree?"] = "Yes" + f" (n = {matching.shape[0]})"
different = hap_phased[hap_phased["likely_denovo_size"] != hap_phased["likely_denovo_size_parsimony"]]
different["Both methods agree?"] = "No" f" (n = {different.shape[0]})"

print (different)


both = pd.concat([different, matching])

both["is_exp"] = both["likely_denovo_size_parsimony"] > 0
conting = both.groupby(["Both methods agree?", "is_exp"]).size().reset_index().rename(columns={0: "count"})
print (conting)
print (conting, conting["count"], conting["count"].values.reshape((2, 2)))
print (ss.chi2_contingency(conting["count"].values.reshape((2, 2))))

f, ax = plt.subplots(figsize=(8, 6))

sns.scatterplot(
    data=both,
    x="likely_denovo_size",
    y="likely_denovo_size_parsimony",
    hue="Both methods agree?",
    ax=ax,
    palette="colorblind",
    alpha=0.5
)
ax.axline((0, 0), slope=1, ls=":", c="gainsboro")
ax.set_xlabel("DNM size (inferred from parental haplotype)")
ax.set_ylabel("DNM size (inferred by allele length parsimony)")
ax.set_yscale("symlog")
ax.set_xscale("symlog")
sns.despine(ax=ax)
f.tight_layout()
f.savefig(snakemake.output.png, dpi=200)
