from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from collections import defaultdict, Counter
import matplotlib.patches as patches
from scipy.spatial.distance import hamming
from sklearn.neighbors import NearestNeighbors


plt.rc("font", size=18, )

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

TRID = "chr8_2623352_2623487_trsolve"


GENOME = "GRCh38"
COHORT = "hprc"
PALETTE = "colorblind"

hapseq  = []
for record in SeqIO.parse(f"trviz/{COHORT}/{TRID}.{GENOME}_alignment_output.fa", "fasta"):
    hapseq.append(list(record.seq))

hapseq = pd.DataFrame(hapseq)

# identify the unique characters in the encoded
# motif arrays, including "-" gaps
uniq_chars = np.unique(hapseq.values)
print (sorted(uniq_chars))

# map characters to colors
cmap = sns.color_palette(PALETTE, uniq_chars.shape[0])
cdict = dict(zip(range(len(uniq_chars)), cmap))
cdict_int = dict(zip(uniq_chars, range(len(uniq_chars))))


NAN_VAL = -999

cdict_int["-"] = NAN_VAL

n_haps, n_motifs = hapseq.shape

hapseq_asint = hapseq.replace(cdict_int)


mb = NearestNeighbors(n_neighbors=n_haps, metric='euclidean').fit(hapseq_asint)
v = mb.kneighbors(hapseq_asint)
smallest = np.argmin(v[0].sum(axis=1))

idx2 = hapseq.apply(lambda row: row.values[row.values != "-"].shape[0], axis=1)

idx1 = v[1][smallest]

hapseq_asint["similarity"] = idx1
hapseq_asint["length"] = idx2


hapseq_asint = hapseq_asint.sort_values(["length", "similarity"]).drop(columns=["length", "similarity"])

f, ax = plt.subplots(figsize=(8, int(0.1 * n_haps)))

cmap = sns.color_palette("colorblind", len(cdict_int) - 1, as_cmap=True)

hapseq_asint = hapseq_asint[hapseq_asint.columns[::-1]]
hapseq_asint[hapseq_asint == NAN_VAL] = np.nan
sns.heatmap(data=hapseq_asint, ax=ax, cmap=cmap, cbar=False)
for i in range(n_haps):
    ax.axhline(y=i, ls="-", c="gainsboro", lw=0.5)
ax.set_xlabel("Number of motifs on haplotype")
ax.set_ylabel(f"{COHORT.upper()} haplotype (n = {n_haps})")
ax.set_xticks(range(0, n_motifs, 5))
ax.set_xticklabels(range(0, n_motifs, 5))
ax.set_yticks([])
f.tight_layout()
f.savefig(f"{COHORT}.excel.png", dpi=200)
