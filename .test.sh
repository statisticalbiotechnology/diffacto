#!/usr/bin/env bash
cd example
diffacto -i iPRG.novo.pep.csv -samples iPRG.samples.lst -out iPRG.denovo.protein.txt \
  -mc_out iPRG.denovo.protein.FDR -min_samples 4 -impute_threshold 0.9 \
  -use_unique True -log2 False
protfile_size=$(wc -l < iPRG.denovo.protein.txt)
fdrfile_size=$(wc -l < iPRG.denovo.protein.FDR)
if [ $protfile_size -ge '200' ] && [ $fdrfile_size -ge '200' ]; then
    # All OK
    exit 0
else
    # Something is wrong
    exit 1
fi
