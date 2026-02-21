/*
 * =============================================================================
 * Module 5: Alignment Trimming (TrimAl)
 * =============================================================================
 * Removes poorly aligned regions and gap-rich columns.
 * Saves column mapping for position tracking in downstream analyses.
 *
 * Nextflow concepts learned:
 *   - Multiple output files from a single process
 *   - Optional outputs
 * =============================================================================
 */

process TRIM_ALIGNMENT {

    tag "${gene_name}"
    label 'process_low'

    publishDir "${params.outdir}/data/trimmed", mode: 'copy'

    input:
    tuple val(gene_name), path(alignment)

    output:
    tuple val(gene_name), path("${gene_name}.trimmed.fasta"),  emit: alignment
    tuple val(gene_name), path("${gene_name}.colnumbering"),   emit: col_mapping

    script:
    """
    trimal \
        -in ${alignment} \
        -out ${gene_name}.trimmed.fasta \
        -automated1 \
        -colnumbering > ${gene_name}.colnumbering

    echo "[trim] Input columns: \$(head -2 ${alignment} | tail -1 | wc -c)"
    echo "[trim] Output columns: \$(head -2 ${gene_name}.trimmed.fasta | tail -1 | wc -c)"
    """

    stub:
    """
    cp ${alignment} ${gene_name}.trimmed.fasta
    echo "1,2,3" > ${gene_name}.colnumbering
    """
}
