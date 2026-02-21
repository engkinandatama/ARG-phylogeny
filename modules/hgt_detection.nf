/*
 * =============================================================================
 * Module 12: HGT Detection
 * =============================================================================
 * Detects horizontal gene transfer via:
 *   A) Phylogenetic incongruence (gene tree vs species tree)
 *   B) Parametric methods (GC content, Codon Adaptation Index)
 *
 * Nextflow concepts learned:
 *   - Multi-input processes (gene tree + species tree + sequences)
 *   - .join() for matching gene names across channels
 * =============================================================================
 */

process HGT_DETECTION {

    tag "${gene_name}"
    label 'process_medium'

    publishDir "${params.outdir}/results/hgt", mode: 'copy'

    input:
    tuple val(gene_name), path(gene_tree), path(clean_fasta)
    path species_tree

    output:
    tuple val(gene_name), path("${gene_name}_hgt_report.csv"), emit: report
    tuple val(gene_name), path("${gene_name}_hgt_detail.json"), emit: detail

    script:
    """
    tree_comparison.py \
        --gene-tree ${gene_tree} \
        --species-tree ${species_tree} \
        --sequences ${clean_fasta} \
        --gene ${gene_name} \
        --out-csv ${gene_name}_hgt_report.csv \
        --out-json ${gene_name}_hgt_detail.json
    """

    stub:
    """
    echo "gene,rf_distance,au_pvalue,gc_content,cai_score,hgt_verdict" > ${gene_name}_hgt_report.csv
    echo "${gene_name},0.5,0.03,52.1,0.65,Moderate" >> ${gene_name}_hgt_report.csv
    echo '{}' > ${gene_name}_hgt_detail.json
    """
}
