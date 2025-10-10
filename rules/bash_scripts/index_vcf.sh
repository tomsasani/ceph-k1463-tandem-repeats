#!/bin/bash
set -e

module load bcftools

bcftools index --tbi --threads ${snakemake[threads]} ${snakemake_input[vcf]}