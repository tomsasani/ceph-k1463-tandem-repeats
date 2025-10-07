import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import scipy.stats as ss
import tqdm


CI = 95
N_BOOTS = 1_000

plt.rc("font", size=14)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

def get_motif_types(row):
    if row["n_motifs"] > 1:
        if row["max_motiflen"] > 6 and row["min_motiflen"] <= 6:
            return "complex"
        else:
            if row["min_motiflen"] == 1:
                return "homopolymer"
            elif row["min_motiflen"] > 1 and row["max_motiflen"] <= 6:
                return "non-homopolymer STR"
            else:
                return "VNTR"
    else:
        if row["min_motiflen"] == 1:
            return "homopolymer"
        elif 1 < row["min_motiflen"] <= 6:
            return "non-homopolymer STR"
        else:
            return "VNTR"

def check_phase(p):
    if p == "unknown": return p
    else:
        support = float(p.split(":")[1])
        if support < 0.75:
            return "unknown"
        else:
            return p.split(":")[0]


def bootstrap(vals, n: int, bins):
    boots = np.zeros(shape=(n, bins.shape[0] - 1))
    for boot in range(n):
        # random sample
        boot_sizes = np.random.choice(vals, size=vals.shape[0], replace=True)
        hist, edges = np.histogram(boot_sizes, bins=bins)
        hist_fracs = hist / np.sum(hist)
        boots[boot, :] = hist_fracs

    mean = np.mean(boots, axis=0)
    lo_bound = (100 - CI) / 2
    lo = np.percentile(boots, lo_bound, axis=0)
    hi = np.percentile(boots, 100 - lo_bound, axis=0)

    return edges, mean, lo, hi

def get_size_in_motifs(row):
    size = row["likely_denovo_size"]
    motif_size = row["min_motiflen"]

    if size % motif_size != 0:
        return np.nan
    else:
        return size // motif_size

ASSEMBLY = "CHM13v2"


mutations = pd.read_csv("CHM13v2.filtered.tsv", sep="\t", dtype={"paternal_id": str})


# read in akshay decomposition
akshay = pd.read_csv("data/K1463.CHM13v2.DNMs.416.demintr.output", sep="\t")
n_events = akshay.groupby("ID").size().reset_index().rename(columns={0: "n_events"})
akshay = akshay.merge(n_events)

akshay = akshay[akshay["n_events"] == 1]

akshay["akshay_motif"] = akshay["Info"].apply(lambda i: i.split(";")[5].split(":")[1] if len(i.split(";")[5].split(":")) > 1 else "UNK")
akshay["akshay_motif_size"] = akshay["Info"].apply(lambda i: len(i.split(";")[5].split(":")[1]) if len(i.split(";")[5].split(":")) > 1 else "UNK")
akshay = akshay[akshay["akshay_motif"] != "UNK"]

# get mutations with a single motif, these are easy
single_motif_mutations = mutations[mutations["n_motifs"] == 1]

# merge multi-motif mutations with akshays file to decompose
multi_motif_mutations = mutations[mutations["n_motifs"] > 1]
multi_motif_mutations = akshay.merge(multi_motif_mutations, left_on="ID", right_on="trid")

# subset to relevant columns
single_motif_mutations = single_motif_mutations[
    ["likely_denovo_size", "motifs", "min_motiflen", "phase"]
].rename(columns={"motifs": "motif", "min_motiflen": "motif_size"})

multi_motif_mutations = multi_motif_mutations[
    ["likely_denovo_size", "akshay_motif", "akshay_motif_size", "phase"]
].rename(columns={"akshay_motif": "motif", "akshay_motif_size": "motif_size"})

combined = pd.concat([single_motif_mutations, multi_motif_mutations])
# combined = single_motif_mutations

combined["motif_change"] = combined["likely_denovo_size"] // combined["motif_size"]
combined["motif_change_leftover"] = combined["likely_denovo_size"] % combined["motif_size"]

combined = combined[(combined["motif_change_leftover"] == 0) & (combined["motif_change"] != 0)]
print (combined.shape)

# figure out "true" TR type
combined["TR type"] = combined["motif_size"].apply(
    lambda m: (
        "homopolymer" if m == 1 else "non-homopolymer STR" if 2 <= m <= 6 else "VNTR"
    )
)


type2idx = dict(zip(combined["TR type"].unique(), range(combined["TR type"].nunique())))

colors = sns.color_palette("colorblind6", 3)

mutations = mutations[mutations["haplotype_in_parent_consensus"] != "UNK"]

mutations["precursor_AL"] = mutations["precursor_sequence_in_parent"].apply(lambda s: len(s))
mutations["is_exp"] = mutations["likely_denovo_size"] > 0
print (mutations.query("precursor_AL < 10"))
f, ax = plt.subplots()
sns.boxplot(data=mutations.query("precursor_AL >= 10"), x="is_exp", y="precursor_AL", ax=ax)
ax.set_yscale("log")
f.savefig("o.png")