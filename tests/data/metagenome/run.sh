fetchMGs  extraction input/metagenomes metagenome output_multi -t 4
gzip output_multi/*

fetchMGs  extraction input/WANG22-1_SAMN08714533-S013_METAG.scaffolds.min500.sub.fasta.gz metagenome output_WANG22-1_SAMN08714533-S013_METAG -t 4
gzip fetchMGs output_WANG22-1_SAMN08714533-S013_METAG/*


