from collections import defaultdict
import glob


ORIG_PATH = "/scratch/ucgd/lustre-labs/quinlan/data-shared/datasets/Palladium/hifi-bams"

# map samples to list of sub-bams
# samples, seqruns, bcids = [], [], []
SMP2SEQRUNS, SMP2BCIDS = defaultdict(list), defaultdict(list)
for fh in glob.glob(f"{NEW_PATH}/K**/*.bam"):
    sample = fh.split("/")[-2].split("_")[1]
    seqrun = fh.split("/")[-1].split(".")[0]
    bcid = fh.split("/")[-1].split(".")[-2]
    SMP2SEQRUNS[sample].append(seqrun)
    SMP2BCIDS[sample].append(bcid)


def get_orig_sample_bam(wildcards):
    assembly_adj = wildcards.ASSEMBLY.lower().rstrip("v2") if wildcards.ASSEMBLY == "CHM13v2" else wildcards.ASSEMBLY
    return f"{ORIG_PATH}/{wildcards.ASSEMBLY.rstrip("v2")}/{SMP2ALT[wildcards.SAMPLE]}.{assembly_adj}.haplotagged.bam"
    

def get_raw_ccs_bams(wildcards):
    seqruns = SMP2SEQRUNS[wildcards.SAMPLE]
    bcids = SMP2BCIDS[wildcards.SAMPLE]
    o = []
    for s,b in zip(seqruns, bcids):
        o.append(f"data/bam/raw/{wildcards.ASSEMBLY}/{wildcards.SAMPLE}.{s}.{b}.sorted.bam")
    return o


rule align_ccs_bam:
    input:
        ref = lambda wildcards: ASSEMBLY2REF[wildcards.ASSEMBLY],
        bam = NEW_PATH + "/K1463_{SAMPLE}/{SEQRUN}.hifi_reads.{BCID}.bam",
        pbmm2 = "/uufs/chpc.utah.edu/common/HIPAA/u1006375/bin/pbmm2"
    output:
        bam = "data/bam/raw/{ASSEMBLY}/{SAMPLE}.{SEQRUN}.{BCID}.sorted.bam",
        bai = "data/bam/raw/{ASSEMBLY}/{SAMPLE}.{SEQRUN}.{BCID}.sorted.bam.bai",
    threads: 8
    params:
        alt_sample_id = lambda wildcards: SMP2ALT[wildcards.SAMPLE]
    resources:
        mem_mb = 64_000
    script:
        "bash_scripts/align_bam.sh"


rule combine_sample_bams:
    input:
        new_bams = get_raw_ccs_bams,
        old_bam = get_orig_sample_bam
    output:
        bam_fh = NEW_PATH + "/merged/{ASSEMBLY}/{SAMPLE}.merged.bam",
    threads: 8
    shell:
        """
        module load samtools
        
        samtools merge --threads {threads} -O BAM -o {output.bam_fh} {input.new_bams} {input.old_bam}
        """


rule index_sample_bams:
    input:
        bam_fh = NEW_PATH + "/merged/{ASSEMBLY}/{SAMPLE}.merged.bam"
    output:
        bam_idx = NEW_PATH + "/merged/{ASSEMBLY}/{SAMPLE}.merged.bam.bai",
    threads: 4
    shell:
        """
        module load samtools
        
        samtools index -@ {threads} {input.bam_fh}
        """