/*
 * =============================================================================
 * Module 10: Multiple Testing Correction
 * =============================================================================
 * Applies Benjamini-Hochberg (FDR) or Holm-Bonferroni correction to
 * p-values from selection analysis results.
 *
 * Nextflow concepts learned:
 *   - .collect() to gather all results before processing
 *   - .groupTuple() to group by gene
 * =============================================================================
 */

process MULTIPLE_TESTING {

    tag "${gene_name}"
    label 'process_low'

    publishDir "${params.outdir}/results/selection_corrected", mode: 'copy'

    input:
    tuple val(gene_name), path(json_files)

    output:
    tuple val(gene_name), path("${gene_name}_corrected_summary.csv"), emit: corrected
    tuple val(gene_name), path("${gene_name}_significant_sites.csv"), emit: significant

    script:
    def method = params.correction_method ?: 'BH'
    def p_thresh = params.p_threshold ?: 0.05
    """
    correct_pvalues.py \
        --input-dir . \
        --gene ${gene_name} \
        --method ${method} \
        --p-threshold ${p_thresh} \
        --out-summary ${gene_name}_corrected_summary.csv \
        --out-significant ${gene_name}_significant_sites.csv
    """

    stub:
    """
    echo "gene,method,site,p_value,q_value,significant" > ${gene_name}_corrected_summary.csv
    echo "gene,method,site,q_value" > ${gene_name}_significant_sites.csv
    """
}
