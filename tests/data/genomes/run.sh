fetchMGs  extraction input/genomes genome output_multi -t 4 -v
gzip output_multi/*

for i in $(ls input/*gz); do fetchMGs  extraction $i genome output"_"`echo $i | cut -f 2 -d "/" | sed 's/.fa.gz//g'` -t 4 -v; done

gzip output*/*
