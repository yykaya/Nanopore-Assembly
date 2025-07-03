# End-to-End Nanopore Assembly Pipeline

> A  simple pipeline for basecalling, read correction, assembly, and scaffolding of Nanopore sequencing (R9) data. Here, I demonstrate its use with Arabidopsis thaliana as an example genome.

## Overview

This repository provides a modular workflow to process Nanopore FAST5 files through basecalling, FASTQ conversion, error correction, genome assembly, and scaffolding. Designed for ease of use on HPC clusters with conda environments.

Key steps:
1. **Basecalling & FASTQ conversion**: Convert raw FAST5 to FASTQ via `pod5`, `dorado`, and `samtools`.  
2. **Read correction**: Apply `dorado correct` to generate polished FASTA.  
3. **Draft assembly**: Perform draft assembly using Flye (and optionally use Mabs-Flye to merge/polish multiple assemblies).  
4. **Haplotig purging**: Remove redundant contigs with `purge_dups`.  
5. **Hybrid assembly merge**: Merge the top two assemblies by contiguity (N50) using [Quickmerge](https://github.com/mahulchak/quickmerge).  
6. **Scaffolding**: Use RagTag to improve contiguity against a reference.


## Repository Structure

```plaintext
├── bin/
│   └── fq2fa.pipeline.sh    # Main pipeline script

├── fast5                  # Raw FAST5 files input
├── corrected_fastq        # Basecalled and error-corrected FASTQ
├── draft_assembly         # Initial Flye draft assembly
├── polished_assembly      # Assembly after Mabs-Flye polishing (optional)
├── purged_assembly        # Haplotig-purged assembly
└── scaffolded_assembly    # Final scaffolded assembly
```

## Requirements

- Linux-based system (e.g., HPC cluster)  
- Bash (>=4.4)  
- [Miniconda/Conda](https://docs.conda.io/)  
- `pod5` (for FAST5 conversion)  
- `dorado` (basecalling & correction) — [GitHub](https://github.com/nanoporetech/dorado)  
- `samtools`  
- `flye`  
- `Mabs-Flye` (assembly merging & polishing) — [GitHub](https://github.com/shelkmike/Mabs)  
- `purge_dups` — [GitHub](https://github.com/dfguan/purge_dups)
- `quickmerge` (hybrid assembly merge)
- `ragtag`  
- `nucmer` (from MUMmer)  



## Outputs
This will generate intermediate FASTQ, corrected FASTA, draft and polished -and purged assemblies, and a final scaffolded assembly in `${OUTPUT_DIR}`.
* `raw.fastq`: basecalled reads
* `corrected.fasta`: error-corrected reads
* `assembly/`: draft assembly outputs
* `polished-purged/`: purged (haplotig-reduced) assembly
* `scaffolded/`: final scaffolded assembly

