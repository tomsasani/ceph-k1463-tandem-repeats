import pandas as pd
from trviz.decomposer import Decomposer
from trviz.utils import get_sample_and_sequence_from_fasta
from trviz.main import TandemRepeatVizWorker
from trviz.motif_aligner import MotifAligner
from trviz.motif_encoder import MotifEncoder


recurrent = pd.read_csv(snakemake.input.recurrent_dnms, sep="\t")
recurrent = recurrent[recurrent["trid"] == snakemake.wildcards.TRID]

# get motifs
motifs = recurrent["motifs"].unique()[0].split(",")

# get the sample IDs and sequences from the combined FASTA
sample_ids, tr_sequences = get_sample_and_sequence_from_fasta(snakemake.input.fasta)

# manually decompose each sample
tr_decomposer = Decomposer()
decomposed_motifs = []
for sample, seq in zip(sample_ids, tr_sequences):

    decomposed = tr_decomposer.decompose(seq, motifs)
    decomposed_motifs.append(decomposed)

motif_encoder = MotifEncoder()

encoded_motifs = motif_encoder.encode(decomposed_motifs, motif_map_file=snakemake.output.key_tsv)

out_df = pd.DataFrame({"sample_id": sample_ids, "original_sequence": tr_sequences, "encoded_motifs": encoded_motifs,})
out_df.to_csv(snakemake.output.seq_tsv, sep="\t", index=False)

motif_aligner = MotifAligner()
motif_aligner.align(
    sample_ids=sample_ids,
    encoded_vntrs=encoded_motifs,
    vid=f"{snakemake.wildcards.TRID}.{snakemake.wildcards.ASSEMBLY}",
    output_dir=f"trviz/{snakemake.wildcards.COHORT}",
)

tr_visualizer = TandemRepeatVizWorker()
tr_visualizer.generate_trplot(
    "test",
    sample_ids,
    tr_sequences,
    motifs,
    output_name=snakemake.output.png,
    figure_size=(int(0.125 * len(sample_ids)), 12),
)
