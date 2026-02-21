/*
 * =============================================================================
 * Module 4: Codon Alignment (MAFFT)
 * =============================================================================
 * Protein-guided codon alignment:
 *   1. Translate NT → AA
 *   2. Align protein sequences with MAFFT
 *   3. Back-translate to codon-aware nucleotide alignment
 *
 * Nextflow concepts learned:
 *   - Shell vs script blocks
 *   - Conda/Docker directive usage
 * =============================================================================
 */

process ALIGN_CODON {

    tag "${gene_name}"
    label 'process_medium'

    publishDir "${params.outdir}/data/aln", mode: 'copy'

    input:
    tuple val(gene_name), path(fasta)

    output:
    tuple val(gene_name), path("${gene_name}.codon.fasta"), emit: alignment

    script:
    """
    codon_align.py \
        --input ${fasta} \
        --output ${gene_name}.codon.fasta \
        --gene ${gene_name} \
        --threads ${task.cpus}
    """

    stub:
    """
    cp ${fasta} ${gene_name}.codon.fasta
    """
}
