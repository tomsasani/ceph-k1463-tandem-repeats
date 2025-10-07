import pandas as pd

mine = pd.read_csv("csv/recurrent/CHM13v2.recurrent.tsv", sep="\t")
old = pd.read_excel("41586_2025_8922_MOESM13_ESM.xlsx", sheet_name="Supplementary_Table11", skiprows=1)

print (old.merge(mine, how="outer", indicator=True).query("_merge == 'left_only'"))