# snakeblastn - blastn and alignment in snakemake

This is an automated snakemake pipeline [1] for conducting local “blastn” searches on a custom BLAST database. Both the query and the database file (in fasta format) can be defined by the user in the config.yaml. Furthermore, this pipeline automatically extracts the requested sequences into fasta format files using “bedtools getfasta” [2], controls for reverse complements (that can be obtained during BLASTn search) and finally creates alignment files using MUSCLE [3]. 

##Prerequisites

This pipeline can be run on desktop computers, but several requirements can be necessary:

* Linux or MacOS system
* [BLAST&copy;](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
* [Miniconda or Anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
* [Snakemake 5.18.3+](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

##DAG plot

<img src="https://github.com/maxwagn/snakeblastn/dag.svg">


###References

[1] Köster, J., & Rahmann, S. (2012). Snakemake—a scalable bioinformatics workflow engine. Bioinformatics, 28(19), 2520-2522.

[2] Quinlan, A. R., & Hall, I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), 841-842.

[3] Edgar, R. C. (2004). MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic acids research, 32(5), 1792-1797.


