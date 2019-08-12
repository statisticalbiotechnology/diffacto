Diffacto: Differential Factor Analysis for Comparative Shotgun Proteomics
==========================================================================

Requirements
--------------

Anaconda_ Python3.5+

Packages needed:

- numpy 1.10+
- scipy 0.17+
- pandas 0.18+
- networkx 1.10+
- scikit-learn 0.17+
- pyteomics_ 3.3+

.. _Anaconda: https://www.continuum.io/downloads
.. _pyteomics: https://pythonhosted.org/pyteomics

Installation via ``pip``
*************************

::

  pip install numpy scipy pandas networkx scikit-learn pyteomics

Installation via ``conda``
***************************

::

  conda env create -f environment.yml
  source activate diffacto_35


Usage
-----

::

  run_diffacto.py [-h] -i I [-db [DB]] [-samples [SAMPLES]] [-log2 LOG2]
                       [-normalize {average,median,GMM,None}]
                       [-farms_mu FARMS_MU] [-farms_alpha FARMS_ALPHA]
                       [-reference REFERENCE] [-min_samples MIN_SAMPLES]
                       [-use_unique USE_UNIQUE]
                       [-impute_threshold IMPUTE_THRESHOLD]
                       [-cutoff_weight CUTOFF_WEIGHT] [-fast FAST] [-out OUT]
                       [-mc_out MC_OUT]
  optional arguments:
  -h, --help            show this help message and exit
  -i I                  Peptides abundances in CSV format. The first row
                        should contain names for all samples. The first column
                        should contain unique peptide sequences. Missing
                        values should be empty instead of zeros. (default:
                        None)
  -db [DB]              Protein database in FASTA format. If None, the peptide
                        file must have protein ID(s) in the second column.
                        (default: None)
  -samples [SAMPLES]    File of the sample list. One run and its sample group
                        per line, separated by tab. If None, read from peptide
                        file headings, then each run will be summarized as a
                        group. (default: None)
  -log2 LOG2            Input abundances are in log scale (True) or linear
                        scale (False) (default: False)
  -normalize {average,median,GMM,None}
                        Method for sample-wise normalization. (default: None)
  -farms_mu FARMS_MU    Hyperparameter mu (default: 0.1)
  -farms_alpha FARMS_ALPHA
                        Hyperparameter weight of prior probability (default:
                        0.1)
  -reference REFERENCE  Names of reference sample groups (separated by
                        semicolon) (default: average)
  -min_samples MIN_SAMPLES
                        Minimum number of samples peptides needed to be
                        quantified in (default: 1)
  -use_unique USE_UNIQUE
                        Use unique peptides only (default: False)
  -impute_threshold IMPUTE_THRESHOLD
                        Minimum fraction of missing values in the group.
                        Impute missing values if missing fraction is larger
                        than the threshold. (default: 0.99)
  -cutoff_weight CUTOFF_WEIGHT
                        Peptides weighted lower than the cutoff will be
                        excluded (default: 0.5)
  -fast FAST            Allow early termination in EM calculation when noise
                        is sufficiently small. (default: False)
  -out OUT              Path to output file (writing in TSV format).
  -mc_out MC_OUT        Path to MCFDR output (writing in TSV format).
                        (default: None)


Example
-------

- Peptide abundances recorded in log scale. map peptides to the protein database HUMAN.fa, using GMM (Gaussian Mixture Model) for per-sample normalization, read sample groups in the file sampleLables.txt, and output protein quantification result to the file protein.txt. Peptide abundance will be scaled by comparing average abundances of all samples.

::

  python run_diffacto.py -i peptides.csv -log2 True -db HUMAN.fa -normalize GMM -samples sampleLables.txt -out protein.txt


- Peptide abundances recorded in linear scale, using median abundances for per-sample normalization, read sample groups in the file sampleLables.txt, and output protein quantification result to the file protein.txt. Peptide abundance will be scaled by comparing to average abundances of samples labeled as of Sample1 and Sample3 in the sample list. Use peptides unique to the protein and quantified at least in 20 samples. For a given group of sample, if missing values consist more than 70% of the results, impute missing values at half of the minimum non-missing abundance. Apply sequential Monte Carlo permutation tests and estimate MCFDR for differentially expressed proteins.

::

  python run_diffacto.py -i peptides.csv -out protein.txt -normalize median -samples sampleLables.txt -ref Sample1;Sample3  -use_unique True  -min_samples 20  -impute_threshold 0.7 -mc_out protein.MCFDR.txt
