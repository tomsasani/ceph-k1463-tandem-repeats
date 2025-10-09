#!/bin/bash
set -e        

${snakemake_input[pbmm2]} align \
    --num-threads ${snakemake[threads]} \
    --sort \
    --sort-memory 4G \
    --preset CCS \
    --sample ${snakemake_params[alt_sample_id]} \
    --bam-index BAI \
    --unmapped \
    ${snakemake_input[ref]} \
    ${snakemake_input[bam]} \
    ${snakemake_output[bam]}
