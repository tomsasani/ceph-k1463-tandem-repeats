import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from upsetplot import from_memberships, from_contents, plot
from matplotlib_venn import venn2, venn3
from bx.intervals.intersection import Interval, IntervalTree
from collections import defaultdict
import csv

plt.rc("font", size=12)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]


ASSEMBLY = "CHM13v2"
ORTHOGONAL_TECH = "element"
SAMPLE = "2211"

SAMPLES = ["2189", "2216", "2211", "2212", "2298", "2215", "2217", "2187"]

orig, topup = [], []
for fh in glob.glob("csv/filtered_and_merged/*.tsv"):
    kind = fh.split("/")[-1].split(".")[-2]
    sample = fh.split("/")[-1].split(".")[0]
    if sample not in SAMPLES: continue

    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    if sample in ("2189", "2216"):
        df["extra"] = True
    else:
        df["extra"] = False
    if kind == "ORIGINAL":

        orig.append(df)
    else:
        topup.append(df)


orig = pd.concat(orig)
topup = pd.concat(topup)



jaccard = []
for sample in SAMPLES:
    _orig = set(orig[orig["sample_id"] == sample]["trid"].to_list())
    _topup = set(topup[topup["sample_id"] == sample]["trid"].to_list())

    j = len(_orig.intersection(_topup)) / len(_orig.union(_topup))
    r = sum([o in _topup for o in _orig]) / len(_orig)
    print (r)
    print ([o in _orig for o in _topup])
    for metric, m in zip(("jaccard", "recovered"), (j, r)):
        jaccard.append({"sample": sample, "metric": metric, "value": m, "topped up": sample in ("2216", "2189")})
jaccard = pd.DataFrame(jaccard)

f, ax = plt.subplots(figsize=(8, 6))
sns.barplot(data=jaccard, x="sample", y="value", hue="metric", ax=ax)
sns.despine(ax=ax)
ax.set_ylabel("Jaccard similarity between DNMs\nidentified using new and old pipelines")
ax.set_xlabel("Sample ID")
f.savefig("jaccard.png", dpi=200)
