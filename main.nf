#!/usr/bin/env nextflow

/*
 * =============================================================================
 * ARG Phylogenomics Toolkit
 * =============================================================================
 * Config-driven Nextflow DSL2 pipeline for evolutionary analysis of
 * antibiotic resistance genes (ARGs).
 *
 * Usage:
 *   nextflow run main.nf -params-file conf/project.yaml -profile conda
 *   nextflow run main.nf -stub -profile local
 *
 * Author: engkinandatama
 * =============================================================================
 */

nextflow.enable.dsl = 2

// ─── Load YAML config ───────────────────────────────────────────────────────

import org.yaml.snakeyaml.Yaml

def loadConfig(config_path) {
    def yaml = new Yaml()
    return yaml.load(new File(config_path).text)
}

// ─── Help message ────────────────────────────────────────────────────────────

def helpMessage() {
    log.info """
    =============================================================================
     ARG Phylogenomics Toolkit  v${workflow.manifest.version}
    =============================================================================

     Usage:
       nextflow run main.nf -params-file conf/project.yaml -profile conda

     Profiles:
       local       : Run locally with Conda/Mamba
       conda       : Use Conda environment
       docker      : Use Docker container
       singularity : Use Singularity container
       codespaces  : GitHub Codespaces
       gitpod      : Gitpod workspace
       test        : Mini test dataset

     Options:
       --config_file  : Path to project YAML config [default: conf/project.yaml]
       --outdir       : Output directory [default: results]
       --help         : Show this help message

    =============================================================================
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// ─── Import modules ─────────────────────────────────────────────────────────

include { FETCH_SEQUENCES     } from './modules/fetch_sequences'
include { CLUSTER_CDHIT       } from './modules/cluster_cdhit'
include { QC_FILTER           } from './modules/qc_filter'
include { ALIGN_CODON         } from './modules/align_codon'
include { TRIM_ALIGNMENT      } from './modules/trim_alignment'
include { ALIGNMENT_QC        } from './modules/alignment_qc'
include { BUILD_TREE          } from './modules/build_tree'
include { RECOMBINATION_GARD  } from './modules/recombination_gard'

// Selection methods (7 separate processes for parallel execution)
include { SELECTION_SLAC      } from './modules/selection_hyphy'
include { SELECTION_FEL       } from './modules/selection_hyphy'
include { SELECTION_MEME      } from './modules/selection_hyphy'
include { SELECTION_FUBAR     } from './modules/selection_hyphy'
include { SELECTION_BUSTED    } from './modules/selection_hyphy'
include { SELECTION_ABSREL    } from './modules/selection_hyphy'
include { SELECTION_RELAX     } from './modules/selection_hyphy'

include { MULTIPLE_TESTING    } from './modules/multiple_testing'

// Species tree sub-processes
include { FETCH_HOUSEKEEPING              } from './modules/species_tree'
include { ALIGN_HOUSEKEEPING              } from './modules/species_tree'
include { CONCAT_AND_BUILD_SPECIES_TREE   } from './modules/species_tree'

include { HGT_DETECTION       } from './modules/hgt_detection'
include { VARIANT_NAMING      } from './modules/reporting'
include { DOMAIN_ANNOTATION   } from './modules/domain_annotation'
include { DOMAIN_OVERLAY      } from './modules/domain_annotation'
include { METADATA_NETWORK    } from './modules/reporting'
include { TEMPORAL_SIGNAL     } from './modules/reporting'
include { REPORTING           } from './modules/reporting'

// ─── Main workflow ──────────────────────────────────────────────────────────

workflow {
    // Load project config
    def config = loadConfig(params.config_file)

    log.info """
    =============================================================================
     ARG Phylogenomics Toolkit  v${workflow.manifest.version}
    =============================================================================
     Project    : ${config.project.name}
     Organism   : ${config.project.organism}
     Genes      : ${config.genes.targets.collect { it.name }.join(', ')}
     Housekeep  : ${config.genes.housekeeping.collect { it.name }.join(', ')}
     Output     : ${params.outdir}
     Profile    : ${workflow.profile}
    =============================================================================
    """.stripIndent()

    // ─── PHASE 1: Data Acquisition & Cleaning ───────────────────────────

    // Create gene channel from config
    ch_genes = Channel.from(config.genes.targets)
        .map { gene ->
            tuple(gene.name, gene.query, gene.max_records, config.project.email)
        }

    // Module 1: Fetch sequences from NCBI
    FETCH_SEQUENCES(ch_genes)

    // Module 2: CD-HIT clustering (99% identity)
    CLUSTER_CDHIT(FETCH_SEQUENCES.out.fasta)

    // Module 3: QC filter
    QC_FILTER(CLUSTER_CDHIT.out.fasta)

    // ─── PHASE 2: Alignment & QC ───────────────────────────────────────

    // Module 4: Protein-guided codon alignment
    ALIGN_CODON(QC_FILTER.out.fasta)

    // Module 5: Alignment trimming
    TRIM_ALIGNMENT(ALIGN_CODON.out.alignment)

    // Module 6: Alignment QC + saturation test
    ALIGNMENT_QC(TRIM_ALIGNMENT.out.alignment)

    // ─── PHASE 3: Phylogeny & Selection ─────────────────────────────────

    // Module 7: Gene tree building (IQ-TREE2 + rooting)
    BUILD_TREE(TRIM_ALIGNMENT.out.alignment)

    // Module 8: GARD recombination detection
    //   CRITICAL: Must complete BEFORE selection analysis
    RECOMBINATION_GARD(TRIM_ALIGNMENT.out.alignment)

    // Prepare selection input: combine GARD partitions with gene tree
    // Using original alignment if no partitions detected
    ch_selection_input = RECOMBINATION_GARD.out.original_alignment
        .join(BUILD_TREE.out.tree)
        .map { gene_name, alignment, tree ->
            tuple(gene_name, alignment, tree)
        }

    // Module 9: Selection analysis (7 methods in parallel)
    //   ALL use "Bacterial and Plant Plastid" genetic code (Code 11)
    SELECTION_SLAC(ch_selection_input)
    SELECTION_FEL(ch_selection_input)
    SELECTION_MEME(ch_selection_input)
    SELECTION_FUBAR(ch_selection_input)
    SELECTION_BUSTED(ch_selection_input)
    SELECTION_ABSREL(ch_selection_input)
    SELECTION_RELAX(ch_selection_input)

    // Collect all selection results per gene
    ch_all_selection = SELECTION_SLAC.out.result
        .mix(SELECTION_FEL.out.result,
             SELECTION_MEME.out.result,
             SELECTION_FUBAR.out.result,
             SELECTION_BUSTED.out.result,
             SELECTION_ABSREL.out.result,
             SELECTION_RELAX.out.result)
        .groupTuple()
        .map { gene_name, methods, json_files ->
            tuple(gene_name, json_files.flatten())
        }

    // Module 10: Multiple testing correction
    MULTIPLE_TESTING(ch_all_selection)

    // ─── PHASE 4: Comparative Analysis ──────────────────────────────────

    // Module 11: Species tree from housekeeping genes
    ch_housekeeping = Channel.from(config.genes.housekeeping)
        .map { gene ->
            tuple(gene.name, gene.query, gene.max_records, config.project.email)
        }

    FETCH_HOUSEKEEPING(ch_housekeeping)
    ALIGN_HOUSEKEEPING(FETCH_HOUSEKEEPING.out.fasta)

    ch_hk_alns = ALIGN_HOUSEKEEPING.out.alignment
        .map { gene_name, aln -> aln }
        .collect()

    CONCAT_AND_BUILD_SPECIES_TREE(ch_hk_alns)

    // Module 12: HGT detection
    ch_hgt_input = BUILD_TREE.out.tree
        .join(QC_FILTER.out.fasta)

    HGT_DETECTION(ch_hgt_input, CONCAT_AND_BUILD_SPECIES_TREE.out.tree)

    // Module 13: Variant naming (CARD)
    VARIANT_NAMING(QC_FILTER.out.fasta)

    // ─── PHASE 5: Functional & Network ──────────────────────────────────

    // Module 14: Domain annotation (HMMER)
    DOMAIN_ANNOTATION(TRIM_ALIGNMENT.out.alignment)

    // Module 15: Domain-selection overlay
    ch_overlay_input = DOMAIN_ANNOTATION.out.domains
        .join(MULTIPLE_TESTING.out.significant)

    DOMAIN_OVERLAY(ch_overlay_input)

    // Module 16: Metadata network + co-selection
    ch_all_meta = FETCH_SEQUENCES.out.meta
        .map { gene_name, meta -> meta }
        .collect()

    METADATA_NETWORK(ch_all_meta)

    // Module 17: Temporal signal
    ch_temporal_input = BUILD_TREE.out.tree
        .join(FETCH_SEQUENCES.out.meta)

    TEMPORAL_SIGNAL(ch_temporal_input)

    // ─── PHASE 6: Reporting ─────────────────────────────────────────────

    // Module 18: Per-gene reporting
    ch_report_input = QC_FILTER.out.stats
        .join(ALIGNMENT_QC.out.stats)
        .join(ALIGNMENT_QC.out.saturation)
        .join(MULTIPLE_TESTING.out.corrected)
        .join(HGT_DETECTION.out.report)
        .join(DOMAIN_ANNOTATION.out.domains)
        .join(TEMPORAL_SIGNAL.out.result)

    REPORTING(ch_report_input)
}

// ─── On completion ──────────────────────────────────────────────────────────

workflow.onComplete {
    log.info """
    =============================================================================
     Pipeline completed!
     Status   : ${workflow.success ? 'SUCCESS' : 'FAILED'}
     Duration : ${workflow.duration}
     Output   : ${params.outdir}
    =============================================================================
    """.stripIndent()
}
