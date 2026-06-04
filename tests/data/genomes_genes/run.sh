fetchMGs  extraction input/genome_genes_faa gene output_multi -t 4 -v -d input/genome_genes_fna
gzip output_multi/*

for i in $(ls input/*faa.gz); do fetchMGs  extraction $i gene output"_"`echo $i | cut -f 2 -d "/" | sed 's/.faa.gz//g'` -t 4 -v -d `echo $i | sed 's/.faa.gz/.fna.gz/g'`; done

gzip output*/*

