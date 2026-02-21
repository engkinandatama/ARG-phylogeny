/*
 * =============================================================================
 * Module 6: Alignment QC + Saturation Test
 * =============================================================================
 * Reports alignment statistics and tests for substitution saturation.
 * Saturation at 3rd codon positions can mislead selection analysis.
 *
 * Nextflow concepts learned:
 *   - Processes that produce reports/stats (not just data)
 *   - Using Python scripts for complex analysis
 * =============================================================================
 */

process ALIGNMENT_QC {

    tag "${gene_name}"
    label 'process_low'

    publishDir "${params.outdir}/qc/alignment",  mode: 'copy', pattern: '*_aln_stats.csv'
    publishDir "${params.outdir}/qc/saturation", mode: 'copy', pattern: '*_saturation*'

    input:
    tuple val(gene_name), path(alignment)

    output:
    tuple val(gene_name), path("${gene_name}_aln_stats.csv"),       emit: stats
    tuple val(gene_name), path("${gene_name}_saturation.csv"),      emit: saturation
    tuple val(gene_name), path("${gene_name}_saturation.png"),      emit: saturation_plot, optional: true

    script:
    """
    saturation_test.py \
        --input ${alignment} \
        --gene ${gene_name} \
        --out-stats ${gene_name}_aln_stats.csv \
        --out-saturation ${gene_name}_saturation.csv \
        --out-plot ${gene_name}_saturation.png
    """

    stub:
    """
    echo "gene,n_seqs,aln_length,gap_pct,informative_sites" > ${gene_name}_aln_stats.csv
    echo "${gene_name},10,900,5.2,150" >> ${gene_name}_aln_stats.csv
    echo "gene,codon_pos,iss,iss_c,saturated" > ${gene_name}_saturation.csv
    echo "${gene_name},1,0.3,0.7,false" >> ${gene_name}_saturation.csv
    touch ${gene_name}_saturation.png
    """
}
