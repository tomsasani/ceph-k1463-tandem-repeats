import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
import scipy.stats as ss
import statsmodels.api as sm
import statsmodels.formula.api as smf


def get_motif_types(row: pd.Series):
    if row["max_motiflen"] == 1:
        return "homopolymer"
    else:
        if row["simple_motif_size"] == "STR":
            return "non-homopolymer STR"
        else:
            return row["simple_motif_size"]

pd.set_option("display.precision", 8)
plt.rc("font", size=14)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

mutations = pd.read_csv(snakemake.input.mutations, dtype={"sample_id": str, "paternal_id": str}, sep="\t")

mutations["generation"] = mutations["sample_id"].apply(lambda s: "G4A" if s.startswith("2000") else "G4B" if s.startswith("2001") else "G3")

# get sample IDs so we can filter the denominator files
sample_ids = mutations["sample_id"].unique()
# map alternate (NAXXXX) IDs to original (2189) IDs
alt2orig = dict(zip(mutations["alt_sample_id"], mutations["sample_id"]))

metadata = pd.read_csv(snakemake.input.metadata, dtype={"UGRP Lab ID (archive)": str})

mutations = mutations.merge(metadata, left_on="sample_id", right_on="UGRP Lab ID (archive)")

mutations = mutations[mutations["phase"] != "unknown"]

mutations["TR type"] = mutations.apply(lambda row: get_motif_types(row), axis=1)

# alpha_counts = (
#     mutations.groupby(
#         [
#             "sample_id",
#             "PaAge",
#             "phase",
#         ]
#     )
#     .size()
#     .reset_index()
#     .rename(columns={0: "count"})
# )
# alpha_totals = alpha_counts.groupby(["sample_id", "PaAge"]).agg(total=("count", "sum")).reset_index()
# alpha_counts = alpha_counts.merge(alpha_totals)
# alpha_counts = alpha_counts[alpha_counts["phase"] == "dad"]
# alpha_counts["alpha"] = alpha_counts["count"] / alpha_counts["total"]



sample_counts = (
    mutations.groupby(
        [
            "sample_id",
            "PaAge",
            "MaAge",
            "phase",
            "TR type",
        ]
    )
    .size()
    .reset_index()
    .rename(columns={0: "count"})
)
sample_totals = sample_counts.groupby(["sample_id"]).agg(total=("count", "sum")).reset_index()
sample_counts = sample_counts.merge(sample_totals)

sample_counts.rename(columns={"phase": "Parent-of-origin"}, inplace=True)

f, axarr = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(8, 8))
tr2idx = dict(zip(sample_counts["TR type"].unique(), [(0, 0), (0, 1), (1, 0), (1, 1)]))

colors = sns.color_palette("colorblind", 2)
for tr, tr_df in sample_counts.groupby("TR type"):
    i, j = tr2idx[tr]
    ax = axarr[i, j]
    if i == 1:
        ax.set_xlabel("Parental age")
    if j == 0:
        ax.set_ylabel("Number of DNMs")
    ax.set_title(tr)
    for par_i, (par, par_df) in enumerate(tr_df.groupby("Parent-of-origin")):
        age_col = "PaAge" if par == "dad" else "MaAge"
        formula = f"count ~ {age_col}"
        mod = smf.glm(
            formula=formula,
            data=par_df,
            family=sm.families.Poisson(link=sm.families.links.Identity()),
        ).fit()

        preds = par_df[[age_col]]
        preds["const"] = 1.

        predictions = mod.get_prediction(preds[['const', age_col]])
        df_predictions = predictions.summary_frame(alpha=0.05) # 95% confidence interval
        df_predictions[age_col] = par_df[age_col].values
        df_predictions.sort_values(age_col, ascending=True, inplace=True)
        if par == "dad" and tr == "non-homopolymer STR":
            print (mod.summary())

        ax.scatter(
            par_df[age_col],
            par_df["count"],
            c=colors[par_i],
            ec="w",
            s=50,
            lw=1,
            zorder=1,
            label=par if i == 1 and j == 0 else None,
        )
        ax.fill_between(
            df_predictions[age_col],
            df_predictions["mean_ci_lower"],
            df_predictions["mean_ci_upper"],
            color=colors[par_i],
            alpha=0.5,
            ec="none",
            zorder=-1,
        )
        ax.plot(
            df_predictions[age_col],
            df_predictions["mean"],
            color=colors[par_i],
            zorder=0

        )
    sns.despine(ax=ax)
axarr[1, 0].legend(shadow=True, title="Parent-of-origin")
f.tight_layout()
f.savefig(snakemake.output.png, dpi=200)
