/*
 * =============================================================================
 * Module 2: CD-HIT Clustering
 * =============================================================================
 * Clusters sequences at 99% identity to remove near-identical variants.
 * Reduces dataset size while preserving evolutionary signal.
 *
 * Nextflow concepts learned:
 *   - Conda/container directive
 *   - Accessing params in script block
 *   - Additional output files (cluster report)
 * =============================================================================
 */

process CLUSTER_CDHIT {

    tag "${gene_name}"
    label 'process_low'

    publishDir "${params.outdir}/data/clustered", mode: 'copy'

    input:
    tuple val(gene_name), path(fasta)

    output:
    tuple val(gene_name), path("${gene_name}.nr.fasta"),        emit: fasta
    tuple val(gene_name), path("${gene_name}.nr.fasta.clstr"),  emit: clusters

    script:
    def identity = params.clustering_identity ?: 0.99
    def word_size = params.clustering_word_size ?: 10
    """
    cd-hit-est \
        -i ${fasta} \
        -o ${gene_name}.nr.fasta \
        -c ${identity} \
        -n ${word_size} \
        -T ${task.cpus} \
        -M 0 \
        -d 0

    echo "[cd-hit] Input: \$(grep -c '>' ${fasta}) sequences"
    echo "[cd-hit] Output: \$(grep -c '>' ${gene_name}.nr.fasta) representative sequences"
    """

    stub:
    """
    cp ${fasta} ${gene_name}.nr.fasta
    echo "Cluster 0" > ${gene_name}.nr.fasta.clstr
    """
}
