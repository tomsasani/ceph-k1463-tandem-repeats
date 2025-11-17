from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from collections import defaultdict, Counter
import matplotlib.patches as patches
import glob

plt.rc("font", size=18, )

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

TRID = "chr8_2623352_2623487_trsolve"

ASSEMBLY = "GRCh38"
COHORT = "ceph"
PALETTE = "colorblind"

# store unique characters
smps = []

hapseq = []
for record in SeqIO.parse(f"trviz/ceph/{TRID}.{ASSEMBLY}_alignment_output.fa", "fasta"):
    hapseq.append(list(record.seq))
    smps.append(record.id)

hapseq = pd.DataFrame(hapseq)

# figure out how variable each motif across haplotypes
uniq_motifs = np.unique(hapseq)

motif_variability = defaultdict(list)
for i, row in hapseq.iterrows():
    motif_counts = Counter(row.to_list())#.most_common()
    for motif in uniq_motifs:
        if motif == "-": continue
        motif_variability[motif].append(motif_counts[motif])

# remove "static" motifs if desired
static_motifs = []
for motif, cts in motif_variability.items():
    if len(set(cts)) == 1:
        static_motifs.append(motif)

hapseq = hapseq.replace(to_replace={m: "-" for m in static_motifs})

orig_cols = hapseq.columns

# identify the unique characters in the encoded
# motif arrays, including "-" gaps
uniq_chars = np.unique(hapseq.values)

# map characters to colors
cmap = sns.color_palette(PALETTE, uniq_chars.shape[0] - 1)
cdict = dict(zip(uniq_chars[1:], cmap))
cdict[uniq_chars[0]] = "white"

hapseq["sample"] = smps
hapseq["sample_id"] = hapseq["sample"].apply(lambda s: s.split("_")[0])
hapseq["generation"] = hapseq["sample_id"].apply(
    lambda s: (
        "G2A"
        if s == "2209"
        else (
            "G2B"
            if s == "2188"
            else (
                "G4A"
                if (s.startswith("2000") and s != "200080")
                else (
                    "G4B"
                    if (s.startswith("2001") and s != "200100")
                    else "G1" if s in ("2281", "2280", "2213", "2214") else "G3"
                )
            )
        )
    )
)

hapseq["parent_status"] = hapseq["sample_id"].apply(lambda s: (
            "G4A"
            if s in ("2216", "200080") else "G4B" if s in ("2189", "200100")
            else (
                "G3"
                if s in ("2209", "2188")
                else "G2A" if s in ("2280", "2281") else "G2B" if s in ("2214", "2213") else "UNK"
            )
        ))
hapseq["parent_sex"] = hapseq["sample_id"].apply(
    lambda s: (
        "male"
        if s in ["2189", "2214", "2281", "200080", "2209"]
        else "female" if s in ["200100", "2213", "2280", "2216", "2188"] else "UNK"
    )
)

hapseq["has_denovo"] = hapseq["sample"].apply(lambda s: s.split("_")[-1])


# merge hapseq with DNMs that have been annotated with parent and hap of origin
res = []
for fh in glob.glob(f"csv/filtered_and_merged/*.{ASSEMBLY}.TOPUP.tsv"):
    df = pd.read_csv(fh, sep="\t", dtype={"sample_id": str})
    res.append(df)
res = pd.concat(res)
res = res[res["trid"] == TRID][["sample_id", "phase", "haplotype_in_parent_consensus"]]

hapseq = hapseq.merge(res, how="left").fillna({"phase": "UNK", "haplotype_in_parent_consensus": "UNK"})
print (hapseq)

# don't plot G1 as a separate generation
gen_to_plot = [g for g in hapseq["generation"].unique() if g != "G1"]

n_per_gen = hapseq.groupby("generation").size().to_dict()


gs_kw = dict(width_ratios=[1., 1], height_ratios=[1, 1., 3])
f, axarr = plt.subplot_mosaic(
    [
        ["G2A", "G4A"],
        ["G2B", "G4A"],
        ["G3", "G4B"],
    ],
    gridspec_kw=gs_kw,
    figsize=(16, 14),
    layout="constrained",
    sharex=True,
)

def fmt_yticklab(s):
    lab = "    A\n{}\n    B".format(s)

parent_hap2color = {1: "lightblue", 2: "blue", 3: "salmon", 4: "firebrick"}

# loop over the generations
for gen_i, gen in enumerate(["G2A", "G2B", "G3", "G4A", "G4B"]):

    # collect the kids in this generation
    kid_df = hapseq[hapseq["generation"] == gen]

    # ignore spouses outside the family
    kid_df = kid_df[~kid_df["sample_id"].isin(["200080", "200100"])]
    
    # plot the kids second
    kid_df["order"] = 0
    # collect the parents in this generation
    par_df = hapseq[hapseq["parent_status"] == gen]
    # plot the parents first
    par_df["order"] = 1

    gen_df = pd.concat([par_df, kid_df]).sort_values(["order", "parent_sex", "sample_id"])
    # gen_smps = gen_df["sample_id"].to_list()

    order_vals = gen_df["order"].values

    # get a list of haplotypes with DNMs
    has_denovo = gen_df["has_denovo"].values == "denovo"
    gen_df["size"] = gen_df.apply(lambda row: len([m for m in row[orig_cols] if m != "-"]), axis=1)
    sizes = gen_df["size"].to_list()
    # convert characters to integers, including "-",
    # which will be converted to 0s
    gen_hapseq = gen_df[orig_cols].values

    n_haps, n_motifs = gen_hapseq.shape

    # get colors for the specific range of values seen in the hapseq array
    # ax = axarr[gen_i]
    ax = axarr[gen]
    ytick_vals = []

    cur_y = 0
    parent_spaced = False
    for hap_i, hap_seq in enumerate(gen_hapseq):
        # manual hacking to separate parents and individual samples
        # while keeping samples' haplotypes plotted together
        if order_vals[hap_i] == 1 and not parent_spaced:
            cur_y += 1.5
            parent_spaced = True
        if hap_i % 2 == 1:
            cur_y += 0.25
        if hap_i % 2 == 0:
            cur_y += 1
        
        ytick_vals.append(cur_y)

        if hap_i >= gen_hapseq.shape[0] - 4:
            adj_i = gen_hapseq.shape[0] - hap_i
            ax.axhline(cur_y, xmin=0, xmax=(sizes[hap_i] + 7) / n_motifs, c=parent_hap2color[adj_i], lw=6, zorder=0)

        # figure out "end" of allele in MSA
        motifs_seen = False
        for motif_i, motif in enumerate(hap_seq):
            if motif != "-": 
                motifs_seen = True
            if motif == "-" and not motifs_seen: continue
            # print (motif, motifs_seen)
            pos = n_motifs - motif_i
            
            rect = patches.FancyBboxPatch(
                (pos, cur_y - 0.2),
                0.33,
                0.33,
                linewidth=1,
                edgecolor="black",
                facecolor=cdict[motif],
                clip_on=False,
                boxstyle="round",
                zorder=1
            )
            ax.add_patch(rect)
        cur_y += 1
    ax.set_xlim(0, n_motifs + 1.5)
    ax.set_ylim(0, max(ytick_vals) + 1.5)

    ax.set_yticks(ytick_vals)
    yticklabs = gen_df["sample_id"].to_list()
    new_yticklabs = []
    for li, l in enumerate(yticklabs):
        if li % 2 == 0:
            l_fmt = f"{l}-A"
        else:
            l_fmt = f"{l}-B"
        new_yticklabs.append(l_fmt)
    ax.set_yticklabels(new_yticklabs, rotation=0)
    ax.set_title(gen)

    sns.despine(ax=ax)

    for i, t in enumerate(ax.yaxis.get_ticklabels()):
        smp = t.get_text().split("-")[0]
        if has_denovo[i] and smp not in par_df["sample_id"].to_list():
            # figure out which parental haplotype it came from
            parent_info = res[res["sample_id"] == smp]
            phase, hap = parent_info["phase"].values[0], parent_info["haplotype_in_parent_consensus"].values[0].split(":")[0]
            if parent_info["phase"].values[0] == "dad":
                if hap == "B":
                    color = "lightblue"
                elif hap == "A":
                    color = "blue"
            if parent_info["phase"].values[0] == "mom":
                if hap == "B":
                    color = "salmon"
                elif hap == "A":
                    color = "red"
            ax.axhline(y=t.get_position()[1], xmax=(sizes[i] + 7) / n_motifs, lw=6, c=color, zorder=0)
            # t.set_color("red")
            t.set_fontweight("bold")
            #t.set_fontstyle("oblique")
        # set colors of parental haplotypes

        if smp in par_df["sample_id"].to_list():
            t.set_fontweight("bold")
            #if t.get_text().split("-")[0] 


# make custom legend
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

# get colors used in palette
custom_legend = []
custom_legend.append(mpatches.Patch(color="white", ec="k", lw=1))
for char, color in zip(uniq_chars, cmap):
    custom_legend.append(mpatches.Patch(color=color, ec="k", lw=1))

key = pd.read_csv(f"trviz/ceph/{TRID}.{ASSEMBLY}.key.tsv", sep="\t", names=["motif", "encoded_motif", "count"])

key = key[~key["encoded_motif"].isin(static_motifs)]
encoded2motif = dict(zip(key["encoded_motif"], key["motif"]))


axarr["G4B"].legend(
    custom_legend,
    ["other (non-polymorphic)"] + key["motif"].to_list(),
    prop={"family": "monospace", "size": 16},
    frameon=True,
    shadow=True,
    loc="lower right",
    title="Nucleotides in motif"
    # bbox_to_anchor=(0.5, 1.05)
)

axarr["G4B"].set_xlabel("Motif number")
f.tight_layout()
f.savefig("trviz_haplotypes.png", dpi=200)
