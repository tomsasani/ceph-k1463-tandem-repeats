import csv
import gzip

outfh = open(snakemake.output.trgt, "w")

with gzip.open(snakemake.input.gangstr, "rt") as infh:
    csvf = csv.reader(infh, delimiter="\t")
    header = ["chrom", "start", "end", "motif_size", "motif", "structure"]
    for l in csvf:
        d = dict(zip(header, l))
        out = [d["chrom"], d["start"], d["end"]]
        ident = "_".join(out) + "_gangstr"
        info = f"ID={ident};MOTIFS={d['motif']};STRUC=({d['motif']})n"
        out.append(info)
        print ("\t".join(out), file=outfh)
outfh.close()