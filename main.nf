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

// TODO: Import modules as they are built
// include { FETCH_SEQUENCES } from './modules/fetch_sequences'
// include { CLUSTER_CDHIT   } from './modules/cluster_cdhit'
// include { QC_FILTER        } from './modules/qc_filter'
// include { ALIGN_CODON      } from './modules/align_codon'
// include { TRIM_ALIGNMENT   } from './modules/trim_alignment'
// include { ALIGNMENT_QC     } from './modules/alignment_qc'
// include { BUILD_TREE       } from './modules/build_tree'
// include { RECOMBINATION_GARD } from './modules/recombination_gard'
// include { SELECTION_HYPHY  } from './modules/selection_hyphy'
// include { MULTIPLE_TESTING } from './modules/multiple_testing'
// include { SPECIES_TREE     } from './modules/species_tree'
// include { HGT_DETECTION    } from './modules/hgt_detection'
// include { VARIANT_NAMING   } from './modules/variant_naming'
// include { DOMAIN_ANNOTATION } from './modules/domain_annotation'
// include { DOMAIN_OVERLAY   } from './modules/domain_overlay'
// include { METADATA_NETWORK } from './modules/metadata_network'
// include { TEMPORAL_SIGNAL  } from './modules/temporal_signal'
// include { REPORTING        } from './modules/reporting'

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
     Output     : ${params.outdir}
     Profile    : ${workflow.profile}
    =============================================================================
    """.stripIndent()

    // Create gene channel from config
    // ch_genes = Channel.from(config.genes.targets)
    //     .map { gene -> tuple(gene.name, gene.query, gene.max_records) }

    // TODO: Wire modules as they are built
    // Pipeline flow:
    //   1. FETCH_SEQUENCES(ch_genes)
    //   2. CLUSTER_CDHIT(FETCH_SEQUENCES.out.fasta)
    //   3. QC_FILTER(CLUSTER_CDHIT.out.fasta)
    //   4. ALIGN_CODON(QC_FILTER.out.fasta)
    //   5. TRIM_ALIGNMENT(ALIGN_CODON.out.alignment)
    //   6. ALIGNMENT_QC(TRIM_ALIGNMENT.out.alignment)
    //   7. BUILD_TREE(TRIM_ALIGNMENT.out.alignment)
    //   8. RECOMBINATION_GARD(TRIM_ALIGNMENT.out.alignment)
    //   9. SELECTION_HYPHY(RECOMBINATION_GARD.out.partitions, BUILD_TREE.out.tree)
    //  10. MULTIPLE_TESTING(SELECTION_HYPHY.out.results)
    //  11. SPECIES_TREE(config.genes.housekeeping)
    //  12. HGT_DETECTION(BUILD_TREE.out.tree, SPECIES_TREE.out.tree)
    //  13. VARIANT_NAMING(QC_FILTER.out.fasta)
    //  14. DOMAIN_ANNOTATION(TRIM_ALIGNMENT.out.alignment)
    //  15. DOMAIN_OVERLAY(DOMAIN_ANNOTATION.out.domains, MULTIPLE_TESTING.out.corrected)
    //  16. METADATA_NETWORK(FETCH_SEQUENCES.out.meta)
    //  17. TEMPORAL_SIGNAL(BUILD_TREE.out.tree, FETCH_SEQUENCES.out.meta)
    //  18. REPORTING(all results collected)

    log.info "Pipeline structure initialized. Modules will be added incrementally."
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
