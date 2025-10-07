import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from typing import List
import glob

plt.rc("font", size=10)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]

ASSEMBLY = "CHM13v2"

mutations = pd.read_csv(f"{ASSEMBLY}.filtered.tsv", dtype={"sample_id": str, "paternal_id": str}, sep="\t")

mutations["generation"] = mutations["sample_id"].apply(lambda s: "G4A" if s.startswith("2000") else "G4B" if s.startswith("2001") else "G3")

col = "child_ratio"

f, ax = plt.subplots()
sns.stripplot(data=mutations, x="generation", y=col, hue="generation", ax=ax, alpha=0.5)
sns.boxplot(data=mutations, x="generation", y=col, hue="generation", ax=ax, fliersize=0)
f.savefig("misc.png", dpi=200)
