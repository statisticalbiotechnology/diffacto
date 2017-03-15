## Diffacto: Examples
----


#### Print usage information
<code>python run_diffacto.py  -h </code>


---
#### Example-1:

<code>
python run_diffacto.py -i iPRG.novo.pep.csv -samples iPRG.samples.lst -out iPRG.denovo.protein.txt -mc_out iPRG.denovo.protein.FDR -min_samples 4 -impute_threshold 0.9 -use_unique True -log2 False
</code>

# 

* input-1, peptide abundances: _iPRG.novo.pep.csv_
* input-2, sample list: _iPRG.samples.lst_
* output-1, protein quantification: _iPRG.denovo.protein.txt_
* output-2, FDR estimation by MC tests: _iPRG.denovo.protein.FDR_
* other parameters:  
    -min_samples 4 (peptide quantified in at least four runs)   
    -impute_threshold 0.9 (threshold for missing value imputation 90%)  
    -use_unique True (only use unique peptides for quantification)  
    -log2 False (input abundances are not in log scale)


---
#### Example-2:

<code>
python run_diffacto.py -i HBY20Mix.peptides.csv -log2 False -samples HBY20Mix.samples.lst -db UP000002311_559292.fasta -out HBY20Mix.protein.txt -min_samples 30 -impute_threshold 0.7 -log2 False -reference REF
</code>

#

* input-1, peptide abundances: _HBY20Mix.peptides.csv_
* input-2, sample list: _HBY20Mix.samples.lst_
* input-3, protein database: _UP000002311_559292.fasta_
* output-1, protein quantification: _HBY20Mix.protein.txt_
* other parameters:  
    -min_samples 30 (peptide quantified in at least 30 runs)   
    -impute_threshold 0.7 (threshold for missing value imputation 70%)  
    -log2 False (input abundances are not in log scale)  
    -reference REF (use the runs labeled 'REF' as the internal reference)  
