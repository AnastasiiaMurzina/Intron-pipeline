# Intron-pipeline
Pipeline to get no overlap introns from fasta with gtf annotation

Run exon_counter.py: $python annotation.gtf [feature=exon gene_type=protein_coding] [output_name]

Run seq_distinguisher.py: $python seq_distinguisher.py sequence.fa intron_nooverlaped_direct.tsv intron_nooverlaped_reverse.tsv result.fasta

