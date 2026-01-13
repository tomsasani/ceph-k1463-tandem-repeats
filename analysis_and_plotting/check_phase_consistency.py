import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def check_phase(p, thresh: float = 0.9):
    if "unknown" in p: return "unknown"
    else:
        n_inf, support = list(map(float, p.split(":")[1:]))
        if support < thresh:
            return "unknown"
        else:
            return "maternal" if p.split(":")[0] == "mom" else "paternal"

a = pd.read_csv("csv/phasing/2189.CHM13v2.TOPUP.phased.3gen.tsv", sep="\t")
b = pd.read_csv("csv/phasing/2216.CHM13v2.TOPUP.phased.3gen.tsv", sep="\t")
gen2 = pd.concat([a, b])
both = gen2[(~gen2["phase_consensus"].str.contains("unknown")) & (~gen2["parent_of_origin_3gen"].str.contains("unknown"))]
print (both[["phase_consensus", "parent_of_origin_3gen"]])

res = []
for thresh in np.arange(0.55, 1.05, 0.05):
    thresh = round(thresh, 3)
    for bs in range(100):
        _both = both.sample(frac=1, replace=True)
        
        _both["p2"] = _both["phase_consensus"].apply(lambda p: check_phase(p, thresh=thresh))
        _both["p3"] = _both["parent_of_origin_3gen"].apply(lambda p: p.split(":")[0])
        if bs == 0 and thresh == 1:
            print (_both[["p2", "p3"]])
        non_unk = _both[_both["p2"] != "unknown"]
        # measure consistency assuming 3gen is best
        res.append(
            {
                "thresh": round(thresh, 2),
                "total": non_unk.shape[0],
                "agree": np.sum(non_unk["p2"] == non_unk["p3"]),
                "trial": bs,
            }
        )

res = pd.DataFrame(res)
res["frac"] = res["agree"] / res["total"]
print (res)

f, (ax1, ax2) = plt.subplots(2)
sns.barplot(data=res, x="thresh", y="frac", ax=ax1, errorbar=("ci", 95))
sns.barplot(data=res, x="thresh", y="total", ax=ax2, errorbar=("ci", 95))
f.savefig("o.png")
