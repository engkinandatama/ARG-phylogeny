/*
 * =============================================================================
 * Module 3: QC / Filter / Dedup
 * =============================================================================
 * Filters sequences based on:
 *   - Minimum nucleotide length
 *   - Maximum ambiguity percentage
 *   - Reading frame check (divisible by 3)
 *   - Internal stop codon removal
 *   - Exact sequence deduplication
 *
 * Nextflow concepts learned:
 *   - path() input/output with dynamic naming
 *   - .map() channel operator (used upstream)
 *   - Script arguments from params
 * =============================================================================
 */

process QC_FILTER {

    tag "${gene_name}"
    label 'process_low'

    publishDir "${params.outdir}/data/clean", mode: 'copy'

    input:
    tuple val(gene_name), path(fasta)

    output:
    tuple val(gene_name), path("${gene_name}.clean.fasta"), emit: fasta
    tuple val(gene_name), path("${gene_name}_qc_stats.csv"), emit: stats

    script:
    def min_len = params.qc_min_len ?: 300
    def max_ambig = params.qc_max_ambig ?: 5
    """
    qc_sequences.py \
        --input ${fasta} \
        --output ${gene_name}.clean.fasta \
        --stats ${gene_name}_qc_stats.csv \
        --gene ${gene_name} \
        --min-len ${min_len} \
        --max-ambig ${max_ambig}
    """

    stub:
    """
    cp ${fasta} ${gene_name}.clean.fasta
    echo "gene,input,kept,removed_short,removed_ambig,removed_frame,removed_stop,removed_dup" > ${gene_name}_qc_stats.csv
    echo "${gene_name},1,1,0,0,0,0,0" >> ${gene_name}_qc_stats.csv
    """
}
