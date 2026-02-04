#!/bin/bash
set -e

mkdir -p ${snakemake_params[output_dir]}
cd ${snakemake_params[output_dir]}

${snakemake_input[trgt_denovo_binary]} trio \
    --reference ${snakemake_params[global_pref]}${snakemake_input[reference]} \
    --bed ${snakemake_params[global_pref]}${snakemake_input[repeats]} \
    --father ${snakemake_params[global_pref]}${snakemake_params[dad_pref]} \
    --mother ${snakemake_params[global_pref]}${snakemake_params[mom_pref]} \
    --child ${snakemake_params[global_pref]}${snakemake_params[kid_pref]} \
    --out ${snakemake_params[global_pref]}${snakemake_output} \
    -@ ${snakemake[threads]} \
    -vv

