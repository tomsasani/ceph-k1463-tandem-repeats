import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as ss

plt.rc("font", size=12)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

res = []
for fh in glob.glob("downsampling/csv/filtered_and_merged/*.tsv"):
    kid_dp, mom_dp, dad_dp = list(map(int, fh.split("/")[-1].split(".")[2:5]))
    df = pd.read_csv(fh, sep="\t")
    df["downsampling"] = f"{kid_dp}:{mom_dp}:{dad_dp}"
    res.append(df)
res = pd.concat(res)

res["min_motif_size"] = res["motif_size"].apply(lambda m: str(m) if int(m) <= 6 else "7+")

res["phase"] = res["phase_consensus"].apply(lambda p: p.split(":")[0] if float(p.split(":")[1]) > 0.75 else "unknown")


BOOTSTRAPS = 100

bootstrapped = []
for ds, ds_df in res[
    res["downsampling"].isin(
        ["50:50:10", "50:50:20", "50:50:30", "50:50:40", "50:50:50"]
    )
].groupby("downsampling"):
    for trial in range(BOOTSTRAPS):
        ds_df_resampled = ds_df.sample(frac=1, replace=True)
        ds_df_resampled["t"] = trial
        bootstrapped.append(ds_df_resampled)
bootstrapped = pd.concat(bootstrapped)

grouped = bootstrapped.groupby(["min_motif_size", "t", "downsampling"]).size().reset_index().rename(columns={0: "count"})
totals = grouped.groupby(["t", "downsampling"]).agg(total=("count", "sum")).reset_index()
grouped = grouped.merge(totals)
grouped["frac"] = grouped["count"] / grouped["total"]

f, ax = plt.subplots(figsize=(8, 6))
sns.barplot(data=grouped, x="min_motif_size", y="frac", hue="downsampling", ax=ax)
ax.legend(title="Downsampled depth (kid:mom:dad)", prop={'family':"monospace", "size": 12})
sns.despine(ax=ax)
ax.set_xlabel("Minimum motif size within TR locus")
ax.set_ylabel("Fraction of DNMs")
f.savefig("downsampled.png", dpi=200)


res = res[res["phase"] != "unknown"]
def how_many_downsampled(d):
    kid, mom, dad = list(map(int, d.split(":")))
    return (mom < 50) + (dad < 50)

def relative_downsampling(row):
    kid, mom, dad = list(map(int, row["downsampling"].split(":")))
    if mom < dad:
        return mom / dad
    elif dad < mom:
        return dad / mom
    else: 
        return 1.

res["n_downsampled"] = res["downsampling"].apply(lambda d: how_many_downsampled(d))

# res = res.query("n_downsampled == 1")
res["downsampled_parent"] = res["downsampling"].apply(lambda d: "mom" if int(d.split(":")[1]) < 50 else "dad" if int(d.split(":")[2]) < 50 else "neither")
res["Downsampling relative to other parent"] = res.apply(lambda row: relative_downsampling(row), axis=1)

neither = res[res["downsampled_parent"] == "neither"]
exp_frac = neither.groupby("phase").size()[0] / neither.shape[0]

grouped = (
    res[res["downsampled_parent"] != "neither"].groupby(
        ["phase", "downsampled_parent", "Downsampling relative to other parent"]
    )
    .size()
    .reset_index()
    .rename(columns={0: "count"})
)

totals = grouped.groupby(["downsampled_parent", "Downsampling relative to other parent"]).agg(total=("count", "sum")).reset_index()
grouped = grouped.merge(totals)
grouped["frac"] = grouped["count"] / grouped["total"]
grouped = grouped[grouped["phase"] == "dad"]

f, ax = plt.subplots(figsize=(8, 6))
sns.barplot(data=grouped, x="downsampled_parent", y="frac", hue="Downsampling relative to other parent", ax=ax)
ax.axhline(y=exp_frac, ls=":", c="darkgrey", zorder=-1)
sns.despine(ax=ax)
ax.set_xlabel("Which parent was downsampled?")
ax.set_ylabel("Fraction of DNMs assigned\nto the paternal germline")
f.savefig("downsampled_phase.png", dpi=200)


res = res[res["downsampling"].isin(["50:50:10", "50:50:20", "50:50:30", "50:50:40", "50:50:50"])]

f, ax = plt.subplots(figsize=(8, 6))
res_bs = []
for bs in range(BOOTSTRAPS):
    res_resampled = res.sample(frac=1, replace=True)
    for (ds, phase), ds_df in res_resampled.groupby(["downsampling", "phase"]):
        if phase == "unknown": continue
        count = ds_df.shape[0]
        res_bs.append({"ds": ds, "phase": phase, "count": count})

res_bs = pd.DataFrame(res_bs)

sns.barplot(data=res_bs, x="phase", y="count", hue="ds", ax=ax)

ax.set_ylabel("Number of DNMs")
ax.set_xlabel("Parent-of-origin")
ax.set_title("")
sns.despine(ax=ax)
ax.legend(title="Downsampled depth (kid:mom:dad)", prop={'family':"monospace", "size": 10})
f.savefig("downsampled_sizes.png", dpi=200)
