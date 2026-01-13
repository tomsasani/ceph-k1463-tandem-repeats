from cyvcf2 import VCF
import cyvcf2
import pandas as pd
import tqdm
from typing import List, Dict, Union
from collections import Counter, namedtuple
import numpy as np
from bx.intervals import Interval, IntervalTree

def identify_informative_parent(
    focal: str,
    mom: str,
    dad: str,
    smp2idx: Dict[str, int],
    gts: np.ndarray,
    gq: np.ndarray,
    td: np.ndarray,
    min_gq: int = 20,
    min_dp: int = 10,
):
    """_summary_

    Args:
        focal (peddy.Sample): _description_
        smp2idx (Dict[str, int]): _description_
        gts (np.ndarray): _description_
        gq (np.ndarray): _description_
        td (np.ndarray): _description_
        min_gq (int, optional): _description_. Defaults to 20.
        min_dp (int, optional): _description_. Defaults to 10.

    Returns:
        Union[peddy.Sample, None]: _description_
    """

    mom_idx = smp2idx[mom]
    dad_idx = smp2idx[dad]

    inf_parent = None
    # if both grandparents are available, we can identify informative sites
    # as sites where one grandparent has more ALT alleles than the other. since
    # we require the focal second-generation individual to be HET at any sites,
    # the parent of that focal individual with more ALTs is the one that donated the
    # ALT allele (by definition).

    # if one of the two grandparents has a genotype that doesn't
    # pass muster, skip this informative site.
    bad_parent_gt = False
    for idx in (mom_idx, dad_idx):
        if not gt_ok(idx, gts, gq, td, min_gq=min_gq, min_dp=min_dp):
            bad_parent_gt = True
    if bad_parent_gt:
        return None
    # if both grandparents have good genotypes, the one with more ALTs is the informative one.
    if gts[mom_idx] > gts[dad_idx]:
        inf_parent = mom
    elif gts[dad_idx] > gts[mom_idx]:
        inf_parent = dad

    # if one grandparent is missing, but the other is HOM_REF or HOM_ALT, we can determine
    # which grandparent the second-generational focal individual inherited their ALT from
    if mom_idx is None and dad_idx is not None:
        # if the grandfather has a genotype, but it doesn't pass muster, skip this site
        if gt_ok(dad_idx, gts, gq, td, min_gq=min_gq, min_dp=min_dp):
            if gts[dad_idx] == 2:
                inf_parent = dad
            elif gts[dad_idx] == 0:
                inf_parent = mom
        else:
            return None
    elif dad_idx is None and mom_idx is not None:
        if gt_ok(mom_idx, gts, gq, td, min_gq=min_gq, min_dp=min_dp):
            if gts[mom_idx] == 2:
                inf_parent = mom
            elif gts[mom_idx] == 0:
                inf_parent = dad
        else:
            return None
    
    return inf_parent


# def catalog_gp_ev(
#     gts: np.ndarray,
#     smp2idx: Dict[str, int],
#     mom: peddy.Sample,
#     dad: peddy.Sample,
# ) -> bool:
#     """_summary_

#     Args:
#         gts (np.ndarray): cyvcf2 v.gt_types array
#         smp2idx (Dict[str, int]): dictionary mapping sample IDs to idxs in VCF
#         mom (peddy.Sample): peddy.Sample object for the mother in the trio
#         dad (peddy.Sample): peddy.Sample object for the father in the trio

#     Returns:
#         bool: whether or not there's evidence for the allele of interest in a grandparent
#     """
#     gp_ev = False
#     # loop over the parents of the second-generation parents (i.e.
#     # the first-generation grandparents of the children that share DNMs)
#     for parent in (mom, dad):           
#         # get grandparents of this gen2 individual if available
#         if parent.dad is not None and parent.dad.sample_id in smp2idx:
#             dad_idx = smp2idx[parent.dad.sample_id]
#             if gts[dad_idx] > 0: 
#                 gp_ev = True
#         if parent.mom is not None and parent.mom.sample_id in smp2idx:
#             mom_idx = smp2idx[parent.mom.sample_id]
#             if gts[mom_idx] > 0: 
#                 gp_ev = True    
#     return gp_ev


def variant_pass(
    v: cyvcf2.Variant,
) -> bool:
    """simple utility to make sure a particular variant
    (which we're considering as a candidate informative site)
    is a PASS-ing SNP that isn't in an exclude region

    Args:
        v (cyvcf2.Variant): _description_
        exclude (Dict[IntervalTree]): _description_

    Returns:
        bool: _description_
    """
    if v.var_type != "snp": 
        return False
    if v.FILTER not in ("PASS", None): 
        return False
    
    return True

def ab_ok(idx: int, ab: np.ndarray):
    return 0.2 <= ab[idx] <= 0.8

def gt_ok(
    idx: int,
    gts: np.ndarray,
    gqs: np.ndarray,
    tds: np.ndarray,
    min_gq: int = 20,
    min_dp: int = 10,
) -> bool:
    """simple utility to figure out if a particular
    genotype meets qual and depth filters

    Args:
        idx (int): _description_
        gqs (np.ndarray): _description_
        tds (np.ndarray): _description_
        min_gq (int, optional): _description_. Defaults to 20.
        min_dp (int, optional): _description_. Defaults to 12.

    Returns:
        bool: _description_
    """
    if gts[idx] == 3:
        return False
    # filter on GQ
    if gqs[idx] < min_gq:
        return False
    # filter on depth
    if tds[idx] < min_dp:
        return False
    return True

def catalog_informative_sites(
    *,
    vcf: VCF,
    region: str,
    focal: str,
    dad: str,
    mom: str,
    spouse: str,
    children: List[str],
    kids_with_dnm: List[str],
    smp2idx: Dict[str, int],
    min_gq: int = 20,
    min_dp: int = 10,
    dnm: namedtuple,
):
    """get a list of informative sites with respect to a second-generation individual.
    each site is stored as a namedtuple object with information about its chromosome,
    position, absolute distance (in bp) from the DNM we're trying to phase, as well as the
    sex of the parent/grandparent who transmit the informative alleles.

    informative sites are always collected w/r/t a particular ("focal") second-generation idnividual.

    Args:
        vcf (VCF): _description_
        region (str): _description_
        focal (str): _description_
        spouse (str): _description_
        children (List[str]): _description_
        kids_with_dnm (List[str]): _description_
        smp2idx (Dict[str, int]): _description_
        dnm (namedtuple): _description_
        exclude (Dict[IntervalTree]): _description_
        min_gq (int, optional): _description_. Defaults to 20.
        min_dp (int, optional): _description_. Defaults to 10.

    Returns:
        _type_: _description_
    """

    res = []
    for v in vcf(region):
        # ignore the DNM
        if v.CHROM == dnm.chrom and v.POS - 1 == dnm.start:
            continue
        # only look at SNPs that pass in all members of the family
        if not variant_pass(v):
            continue
        # get metadata for this SNP
        gts = v.gt_types
        ad = v.gt_alt_depths
        td = v.format("DP")[:, 0]
        gq = v.gt_quals
        ab = ad / td

        focal_idx, spouse_idx = smp2idx[focal], smp2idx[spouse]
        # ensure that this is an informative site in the second-generation
        # parents. i.e., one parent is HET and one is HOM_REF.
        if gts[focal_idx] + gts[spouse_idx] != 1:
            continue
        # make sure allele balance is OK in both parents
        bad_gt = False
        for idx in (focal_idx, spouse_idx):
            if not gt_ok(idx, gts, gq, td, min_gq=min_gq, min_dp=min_dp):
                bad_gt = True
        if bad_gt:
            continue

        # make sure the focal parent is HET and the spouse is HOM_REF.
        if not (gts[focal_idx] == 1 and gts[spouse_idx] == 0):
            continue

        inf_parent = identify_informative_parent(
            focal,
            mom,
            dad,
            smp2idx,
            gts,
            gq,
            td,
            min_gq=min_gq,
            min_dp=min_dp,
        )

        # NOTE: do we only *need* informative grandparents when dealing with second-gen DNMs?
        if inf_parent is None:
            continue

        # loop over the third-generation siblings and figure out who
        # inherited the informative allele at this site. we want to be
        # *permissive* when counting children with informative alleles.
        kids_with_informative_allele = []
        for sib in children:
            sib_idx = smp2idx[sib]
            has_inf = (
                gts[sib_idx] == 1
                and gt_ok(
                    sib_idx,
                    gts,
                    gq,
                    td,
                    min_gq=min_gq,
                    min_dp=min_dp,
                )
                and ab_ok(sib_idx, ab)
            )
            if has_inf:
                kids_with_informative_allele.append(sib)

        # we'll only keep track of inheritance patterns in which the informative allele
        # donated by a parent is inherited alongside the DNM (or inherited alone).
        if not all([k in kids_with_informative_allele for k in kids_with_dnm]):
            continue

        # create a list of "statuses" for each child with the informative allele.
        # formatted SAMPLE_ID-Y if they have the DNM or SAMPLE_ID-N if they don't.
        inf_string = []
        for k in kids_with_informative_allele:
            has_dnm = "Y" if k in kids_with_dnm else "N"
            inf_string.append(f"{k}-{has_dnm}")

        if len(inf_string) == 0: 
            continue

        # created named tuple object to store information about this site
        InfSite = namedtuple("InfSite", "chrom pos dist poi focal inh")
        inf_site = InfSite(
            v.CHROM,
            v.POS,
            abs(int(dnm.end) - v.POS),
            inf_parent,
            focal,
            "|".join(inf_string),
        )
        res.append(inf_site)

    return res


def main():
    # assume we have a VCF with SNP genotypes for all individuals in the pedigree
    vcf = VCF(snakemake.input.cohort_vcf, gts012=True)
    smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))

    if snakemake.wildcards.SAMPLE == "2216":
        children = ["200081", "200082", "200084", "200085", "200086", "200087"]
        spouse = "200080"
        focal = "NA12879"

    elif snakemake.wildcards.SAMPLE == "2189":
        children = ["200101", "200102", "200103", "200104", "200105", "200106"]
        spouse = "200100"
        focal = "NA12886"
    
    dad, mom = "NA12877", "NA12878"

    smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))

    mutations = pd.read_csv(snakemake.input.mutation_df, sep="\t")
    # loop over sites that have de novo mutations in at least one individual
    res = []
    for i, row in tqdm.tqdm(mutations.iterrows()):
        row_dict = row.to_dict()

        if type(row["children_with_denovo_allele"]) is float:
            row_dict.update(
                {

                    "parent_of_origin_3gen": "unknown",
                }
            )
            res.append(row_dict)
            continue

        # extract basic information about the DNM and store in a namedtuple object
        DnVar = namedtuple("dnm", "chrom start end")
        chrom = row["#chrom"]
        end = row["end"]
        start = end - 1

        dnm = DnVar(chrom, start, end)

        kids_with_dnm = row["children_with_denovo_allele"].split(",")

        # add some slop to the DNM, around which we'll search for informative sites.
        adj_start, adj_end = end - 1 - snakemake.params.slop, end + snakemake.params.slop
        if adj_start < 0:
            adj_start = 1
        slop_region = f"{chrom}:{adj_start}-{adj_end}"

        par_sites = catalog_informative_sites(
            vcf=vcf,
            region=slop_region,
            focal=focal,
            dad=dad,
            mom=mom,
            spouse=spouse,
            children=children,
            kids_with_dnm=kids_with_dnm,
            smp2idx=smp2idx,
            min_gq=20,
            min_dp=10,
            dnm=dnm,
        )

        if len(par_sites) < 1:
            continue
        else:
            # sort the informative sites by distance to the DNM, and
            # report the 50 closest to the DNM
            sorted_par_sites = sorted(par_sites, key=lambda s: s.dist)[:50]
            # keep track of the "inheritance combos" -- that is, the number of instances where
            # a particular informative allele (known to come from dad or mom) was inherited by
            # a list of children in the family.
            inheritance_combos = [
                ":".join([inf_site.poi, inf_site.focal, inf_site.inh])
                for inf_site in sorted_par_sites
            ]

            most_common_combos = Counter(inheritance_combos).most_common()
            total_combos = len(inheritance_combos)
            

            # calculate the most common combos *for each parent*
            most_common_hap = most_common_combos[0][0]
            most_common_freq = most_common_combos[0][1] / total_combos

            # loop over the kids with the most common combination
            states = []
            for child in children:
                shares_haplotype = False
                shares_dnm = False
                for _sib in most_common_hap.split(":")[-1].split("|"):
                    sibname, has_dnm = _sib.split("-")
                    if sibname == child:
                        shares_haplotype = True
                        if has_dnm == "Y":
                            shares_dnm = True
                    else: continue
                # determine the state of this variant
                state = None
                if shares_haplotype and shares_dnm:
                    state = 2
                elif shares_haplotype and not shares_dnm:
                    state = 1
                elif not shares_haplotype and not shares_dnm:
                    state = 0
                assert not (shares_dnm and not shares_haplotype)
                states.append(state)
            pz = True
            if all([s in (0, 2) for s in states]):
                pz = False
            # if the DNM is in the second generation, the haplotype of origin
            # should be a grandparental haplotype -- if it's in the third generation,
            # the haplotype of origin should be a parental haplotype
            poi_idx = 0
            parent_of_origin = most_common_hap.split(":")[poi_idx]
            poi_sex = "paternal" if parent_of_origin == dad else "maternal"


            row_dict.update(
                {

                    "parent_of_origin_3gen": f"{poi_sex}:{most_common_combos[0][1]}:{most_common_freq}",
                }
            )
            res.append(row_dict)

    res_df = pd.DataFrame(res)

    res_df.to_csv(snakemake.output.tsv, sep="\t", index=False)

if __name__ == "__main__":
    main()
