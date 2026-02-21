/*
 * =============================================================================
 * Module 1: Fetch Sequences from NCBI
 * =============================================================================
 * Downloads CDS sequences from NCBI GenBank using Entrez queries.
 * Outputs raw FASTA + metadata CSV per gene.
 *
 * Nextflow concepts learned:
 *   - process definition (input, output, script)
 *   - val() and path() qualifiers
 *   - publishDir directive
 *   - label for resource allocation
 *   - maxForks to limit concurrent API calls
 * =============================================================================
 */

process FETCH_SEQUENCES {

    tag "${gene_name}"
    label 'process_internet'

    publishDir "${params.outdir}/data/raw",  mode: 'copy', pattern: '*.fasta'
    publishDir "${params.outdir}/data/meta", mode: 'copy', pattern: '*_meta.csv'

    input:
    tuple val(gene_name), val(query), val(max_records), val(email)

    output:
    tuple val(gene_name), path("${gene_name}.fasta"),    emit: fasta
    tuple val(gene_name), path("${gene_name}_meta.csv"), emit: meta

    script:
    """
    fetch_ncbi.py \
        --gene "${gene_name}" \
        --query "${query}" \
        --max-records ${max_records} \
        --email "${email}" \
        --out-fasta "${gene_name}.fasta" \
        --out-meta "${gene_name}_meta.csv"
    """

    stub:
    """
    echo ">stub_seq_1|gene=${gene_name}" > ${gene_name}.fasta
    echo "ATGATGATGATG" >> ${gene_name}.fasta
    echo "acc,gene,product,protein_id,country,date,length_nt" > ${gene_name}_meta.csv
    echo "STUB001,${gene_name},stub_product,STUB_PID,Unknown,2024-01-01,12" >> ${gene_name}_meta.csv
    """
}
