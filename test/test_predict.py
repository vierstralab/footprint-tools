#!/usr/bin/which python3


from genome_tools import genomic_interval, genomic_interval_set

from footprint_tools import cutcounts
from footprint_tools.modeling import predict, bias

from pysam import FastaFile


base_url = "https://resources.altius.org/~jvierstra/projects/footprinting.2020/per.dataset/CD20+-DS18208/"
bam_file_url = base_url + "reads.bam"

bm_file_path = "data/vierstra_et_al.6mer-model.txt"
fasta_file_path = "/home/jvierstra/data/genomes/hg38/hg.all.fa"


interval = genomic_interval('chr19', 48363826, 48364602)

fasta_fn = FastaFile(fasta_file_path)
reads_fn = cutcounts.bamfile(bam_file_url)

bm = bias.kmer_model(bm_file_path)


print(reads_fn[interval])
print(fasta_fn.fetch(interval.chrom, interval.start, interval.end))

pred = predict.prediction(reads_fn, fasta_fn, bm, half_window_width = 5, smoothing_half_window_width = 50, smoothing_clip = 0.01)

obs, exp, win = pred.compute(interval)

print(obs['+'].shape, exp['+'].shape, win['+'].shape)
