# 🧬 ARG Phylogenomics Toolkit

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04-23aa62.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A **config-driven Nextflow DSL2 pipeline** for comprehensive evolutionary analysis of antibiotic resistance genes (ARGs). Designed as a reusable toolkit — configure once via YAML, run for any organism/gene set.

## Overview

This pipeline performs end-to-end phylogenomic analysis of ARGs, from sequence acquisition to publication-ready reports:

```
NCBI Fetch → CD-HIT → QC → Codon Alignment → Trimming → Saturation Test
  → Gene Trees (rooted) → GARD → Selection (HyPhy ×7) → Correction
  → Species Tree → HGT Detection → Domain Overlay → Reporting
```

### Key Features

- **18 analysis steps** in a single automated pipeline
- **Config-driven** — single YAML file controls everything
- **Recombination-aware selection** — GARD partitioning before HyPhy
- **Correct bacterial genetic code** (NCBI Code 11)
- **HGT detection** — phylogenetic incongruence + parametric methods
- **Dual output** — static plots (PNG/PDF) + interactive (HTML)
- **Reproducible** — Conda/Docker containerized
- **Resumable** — Nextflow `-resume` picks up from last checkpoint

## Quick Start

### Prerequisites

- [Nextflow](https://www.nextflow.io/) ≥ 23.04
- [Conda](https://docs.conda.io/) or [Docker](https://www.docker.com/)

### Run

```bash
# 1. Clone
git clone https://github.com/engkinandatama/ARG-phylogeny.git
cd ARG-phylogeny

# 2. Edit config
cp conf/project.yaml conf/my_project.yaml
# Edit genes, organism, parameters...

# 3. Run
nextflow run main.nf -params-file conf/my_project.yaml -profile conda

# 4. Resume if interrupted
nextflow run main.nf -params-file conf/my_project.yaml -profile conda -resume
```

### Test Run

```bash
# Stub run (validates DAG, no data)
nextflow run main.nf -stub -profile local

# Mini test (10 sequences per gene)
nextflow run main.nf -profile test
```

## Pipeline Steps

| # | Step | Tool | Description |
|---|------|------|-------------|
| 1 | Fetch Sequences | BioPython/Entrez | Download CDS from NCBI |
| 2 | Cluster | CD-HIT | Remove near-identical sequences (99%) |
| 3 | QC/Filter | BioPython | Length, ambiguity, frame, stop codons |
| 4 | Codon Alignment | MAFFT | Protein-guided codon alignment |
| 5 | Trim | TrimAl | Remove poorly aligned regions |
| 6 | Alignment QC | Custom | Statistics + Xia's saturation test |
| 7 | Gene Trees | IQ-TREE2 | ML trees with model selection + rooting |
| 8 | Recombination | HyPhy GARD | Breakpoint detection + partitioning |
| 9 | Selection | HyPhy (×7) | SLAC/FEL/MEME/FUBAR/BUSTED/aBSREL/RELAX |
| 10 | Correction | Python | Multiple testing correction (FDR) |
| 11 | Species Tree | MAFFT + IQ-TREE2 | Housekeeping gene concatenation tree |
| 12 | HGT Detection | ete3 + IQ-TREE2 | Topology comparison + GC/CAI analysis |
| 13 | Variant Naming | BLAST + CARD | Standardized ARG variant names |
| 14 | Domain Annotation | HMMER | Pfam domain mapping |
| 15 | Domain Overlay | matplotlib/Plotly | Selection sites on domain architecture |
| 16 | Network | NetworkX | Strain-gene co-occurrence + co-selection |
| 17 | Temporal Signal | Custom | Root-to-tip regression |
| 18 | Reporting | Jinja2 | Per-gene reports + summary + methods |

## Configuration

All pipeline behavior is controlled by a single YAML file. See [`conf/project.yaml`](conf/project.yaml) for the full template with documentation.

```yaml
project:
  name: "ARG-Klebsiella"
  organism: "Klebsiella pneumoniae"

genes:
  targets:
    - name: blaCTX-M
      query: "(Klebsiella pneumoniae[Organism]) AND blaCTX-M[Gene]..."
```

## Current Target

**Antibiotic Resistance Genes in *Klebsiella pneumoniae***:

- Beta-lactamases: blaCTX-M, blaTEM, blaSHV
- Quinolone resistance: oqxAB
- Sulfonamide resistance: sul1
- Tetracycline resistance: tetA

## Output

```
results/
├── qc/                    # QC statistics per gene
├── trees/                 # Gene trees + species tree (Newick)
├── selection/             # HyPhy JSON results (corrected p-values)
├── recombination/         # GARD breakpoints
├── hgt/                   # HGT evidence reports
├── domains/               # Pfam annotations + overlay plots
├── network/               # Co-occurrence graph (GEXF)
├── temporal/              # Root-to-tip regression
├── figures/               # All plots (PNG + HTML)
└── reports/               # Per-gene + summary reports
```

## License

MIT

## Citation

If you use this pipeline in your research, please cite the individual tools used (listed in the auto-generated Methods paragraph in the report output).
