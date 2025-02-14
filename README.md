# ZFY evolutionary analysis

This repo provides analysis scripts for Garcia et al.

## Dependencies

Some software was already installed on our cluster and on the `PATH`, other software was added locally. The scripts expect the following:

### Scripts and binaries

The following binaries must be supplied in `./bin`:

Program             | Source  | Reference  | Version used in this study
--------------------|---------|------------|---------------------------
`macse_v2.07.jar` | https://www.agap-ge2pop.org/macse/ | [doi:10.1371/journal.pone.0022594](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0022594) | 2.07
`muscle5.1.win64.exe` or `muscle5.1.linux_intelx64` | https://drive5.com/muscle5/ | [doi:10.1038/s41467-022-34630-w](https://www.nature.com/articles/s41467-022-34630-w) | 5.1
`nlstradamus.pl` | http://www.moseslab.csb.utoronto.ca/NLStradamus/ | [doi:10.1186/1471-2105-10-202](https://doi.org/10.1186/1471-2105-10-202) | -
`pwm_predict/pwm_predict` | https://zf.princeton.edu/ | [doi:10.1093/nar/gkt890](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3874201/) | -
`geneconv.exe` | https://www.math.wustl.edu/~sawyer/geneconv/ | [10.1093/oxfordjournals.molbev.a040567](https://pubmed.ncbi.nlm.nih.gov/2677599/) | 1.81

Binaries expected on the `PATH`:

Program             | Source  | Reference  | Version used in this study
--------------------|---------|------------|---------------------------
`PAML ` | http://abacus.gene.ucl.ac.uk/software/paml.html | [doi:10.1093/molbev/msm088](https://pubmed.ncbi.nlm.nih.gov/17483113/) | 4.9h
`IQ-TREE` | www.iqtree.org/ | [doi:10.1093/molbev/msu300](https://doi.org/10.1093/molbev/msu300) | 1.6.12
`hmmsearch` | http://hmmer.org/ | [doi:10.1371/journal.pcbi.1002195](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002195) | HMMER v3.3

Other software:

Program             | Source  | Reference  | Version used in this study
--------------------|---------|------------|---------------------------
HyPhy | http://www.hyphy.org/ | [doi:10.1093/bioinformatics/bti079](https://doi.org/10.1093/bioinformatics/bti079) | 2.5.29
Divvier | https://github.com/simonwhelan/Divvier | doi: 10.1093/molbev/msz142(https://pmc.ncbi.nlm.nih.gov/articles/PMC6933875/) | 1.01

HyPhy and Divvier are installed as python packages in conda environments:
```
conda create -n hyphy           # create conda environment
conda activate hyphy            # enter the environment
conda install -c bioconda hyphy # install hyphy from bioconda

conda create -n divvier           # create conda environment
conda activate divvier            # enter the environment
conda install -c bioconda divvier # install divvier from bioconda

# To activate a conda environment from within a shell script:
source activate hyphy

```

### R packages

R packages used are listed in the `src/functions.R` script. If any packages are missing, the script will try to install them.

## Running the analyses

FASTA sequences are provided in the `fasta` directory.


### Generate data files

Invoke the primary analysis via:
```
Rscript alignNucleotide.R
```

This will generate output data in the directory `aln`.

`codeml` takes a long time to run, so `alignNucleotide.R` generates a batch submission script for the cluster instead. Modify this as needed for your own system.

```
bash run_paml.sh
```

### Testing for gene conversion

GENECONV was used to test gene conversion, but it is Windows only (rather, the Unix code did not compile on our cluster), so the code is provided as a separate script.

```
Rscript runGeneconv.R
```

### Create figures

Create the figures from the data in `aln`:

```
Rscript createFigures.R
```

This will output figures to `./figure`

