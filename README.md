# ZFY evolutionary analysis

This repo provides analysis scripts for Garcia et al.

FASTA sequences are provided in the `fasta` directory.

Invoke the primary analysis via:
```
Rscript alignNucleotide.R
```

This will generate output files in the directory `aln`.

`codeml` takes a long time to run, so we generate a batch submission script for the cluster instead. Modify this as needed for your own system

```
bash run_paml.sh
```

Create the figures from the data in `aln`:

```
Rscript createFigures.R
```

Run GENECONV for testing gene conversion (Windows only, so provided as a separate script)

```
Rscript runGeneconv.R
```