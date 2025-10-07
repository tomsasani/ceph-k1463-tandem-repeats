import pandas as pd
import numpy as np
import argparse
from collections import Counter
from typing import List
import matplotlib.pyplot as plt
from schema import HaplotypedSchema
import utils


def measure_consistency(
    df: pd.DataFrame,
    columns: List[str],
    plot: bool = False,
) -> pd.DataFrame:
    """given a dataframe that contains information about informative sites
    surrounding a candidate DNM, find the longest stretch of 'consistent'
    informative sites that all support the same haplotype assignment.

    Args:
        df (pd.DataFrame): pandas DataFrame object containing information about informative sites.
        columns (List[str]): list of columns to use in consistency checks

    Returns:
        pd.DataFrame: pandas DataFrame containing a subset of only the informative sites in the
            longest continuous stretch of consistent sites.
    """

    # sort the informative sites by absolute distance to the STR
    df_sorted = df.sort_values(
        "abs_diff_to_str",
        ascending=True,
    ).reset_index(drop=True)

    sorted_values = df_sorted[columns].values

    # NOTE: this is a new approach. if every one of these sites is in the same
    # phase block, why not just take a consensus approach?
    freqs = Counter(sorted_values).most_common()
    most_common_poi = freqs[0][0]
    poi_frac = freqs[0][1] / len(sorted_values)

    return most_common_poi, len(sorted_values), poi_frac


def main(args):

    dtypes = utils.DTYPES.copy()
    dtypes.update({"sample_id": "string"})

    inf_sites = pd.read_csv(args.annotated_dnms, sep="\t", dtype=dtypes)
    HaplotypedSchema.validate(inf_sites)

    COLS = ["trid", "sample_id", "genotype", "index", "suffix"]

    res = []

    plot = False

    # get a dataframe with N rows for each DNM. each of the N rows
    # represents a single informative site
    for (
        trid,
        sample_id,
        genotype,
        index,
        suffix,
    ), df in inf_sites.groupby(COLS):
                
        # if we're on a male sex chromosome, the phase is simple
        is_male_x = suffix.startswith("S") and trid.startswith("chrX")
        is_male_y = suffix.startswith("S") and trid.startswith("chrY")

        poi, poi_support = "unknown", 0
        hoi, hoi_support = "unknown", 0

        if is_male_x:
            poi = "mom"
            poi_support = 1
        elif is_male_y:
            poi = "dad"
            poi_support = 1
        elif df.shape[0] == 0:
            poi = "unknown"
            poi_support = 0

        # otherwise, we have to use informative sites
        else:

            poi, poi_inf, poi_support = measure_consistency(df, "str_parent_of_origin", plot=True if plot is False else False)
            plot = True
            # these sites have to be the same as the ones for which the parent of origin is consistent
            hap_sites = df[(df["haplotype_in_parent"] != "unknown") & (df["str_parent_of_origin"] == poi)]
            if hap_sites.shape[0] == 0:
                hoi, hoi_inf = "unknown", 0
            else:
                hoi, hoi_inf, hoi_support = measure_consistency(hap_sites, "haplotype_in_parent")

        
        df["phase_consensus"] = f"{poi}:{poi_inf}:{poi_support}"
        df["haplotype_in_parent_consensus"] = f"{hoi}:{hoi_inf}:{hoi_support}"

        res.append(df)

    res_df = (
        pd.concat(res)
        .drop_duplicates(COLS)
        .drop(
            columns=[
                "abs_diff_to_str",
                "inf_chrom",
                "inf_pos",
                "dad_inf_gt",
                "mom_inf_gt",
                "str_parent_of_origin",
                "haplotype_in_parent",
            ]
        )
    )
    print (res_df.shape)

    res_df.to_csv(args.out, sep="\t", index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--annotated_dnms")
    p.add_argument("--out")
    args = p.parse_args()
    main(args)
