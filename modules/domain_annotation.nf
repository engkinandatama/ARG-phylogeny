/*
 * =============================================================================
 * Module 14: Domain Annotation (HMMER)
 * =============================================================================
 */

process DOMAIN_ANNOTATION {
    tag "${gene_name}"
    label 'process_medium'
    publishDir "${params.outdir}/results/domains", mode: 'copy'

    input:
    tuple val(gene_name), path(alignment)

    output:
    tuple val(gene_name), path("${gene_name}_domains.json"), emit: domains

    script:
    """
    domain_overlay.py annotate \
        --input ${alignment} \
        --gene ${gene_name} \
        --output ${gene_name}_domains.json
    """

    stub:
    """
    echo '{"domains": []}' > ${gene_name}_domains.json
    """
}

/*
 * =============================================================================
 * Module 15: Domain-Selection Overlay
 * =============================================================================
 */

process DOMAIN_OVERLAY {
    tag "${gene_name}"
    label 'process_low'
    publishDir "${params.outdir}/figures", mode: 'copy'

    input:
    tuple val(gene_name), path(domains_json), path(selection_csv)

    output:
    tuple val(gene_name), path("${gene_name}_domain_overlay.png"), emit: static_plot
    tuple val(gene_name), path("${gene_name}_domain_overlay.html"), emit: interactive_plot

    script:
    """
    domain_overlay.py plot \
        --domains ${domains_json} \
        --selection ${selection_csv} \
        --gene ${gene_name} \
        --out-png ${gene_name}_domain_overlay.png \
        --out-html ${gene_name}_domain_overlay.html
    """

    stub:
    """
    touch ${gene_name}_domain_overlay.png
    echo "<html>stub</html>" > ${gene_name}_domain_overlay.html
    """
}
