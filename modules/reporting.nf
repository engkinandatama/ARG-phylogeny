/*
 * =============================================================================
 * Module 16: Metadata Network + Co-selection
 * =============================================================================
 */

process METADATA_NETWORK {
    tag "network"
    label 'process_low'
    publishDir "${params.outdir}/results/network", mode: 'copy'

    input:
    path meta_files  // All metadata CSVs collected

    output:
    path "arg_strain_network.gexf",    emit: network
    path "co_selection_matrix.csv",    emit: co_selection
    path "co_selection_heatmap.png",   emit: heatmap, optional: true

    script:
    """
    build_network.py \
        --input-dir . \
        --out-gexf arg_strain_network.gexf \
        --out-coselection co_selection_matrix.csv \
        --out-heatmap co_selection_heatmap.png
    """

    stub:
    """
    echo "stub" > arg_strain_network.gexf
    echo "gene_a,gene_b,jaccard,p_value" > co_selection_matrix.csv
    touch co_selection_heatmap.png
    """
}

/*
 * =============================================================================
 * Module 17: Temporal Signal
 * =============================================================================
 */

process TEMPORAL_SIGNAL {
    tag "${gene_name}"
    label 'process_low'
    publishDir "${params.outdir}/results/temporal", mode: 'copy'

    input:
    tuple val(gene_name), path(tree), path(meta)

    output:
    tuple val(gene_name), path("${gene_name}_temporal.csv"),  emit: result
    tuple val(gene_name), path("${gene_name}_temporal.png"),  emit: plot, optional: true

    script:
    """
    temporal_signal.py \
        --tree ${tree} \
        --meta ${meta} \
        --gene ${gene_name} \
        --out-csv ${gene_name}_temporal.csv \
        --out-plot ${gene_name}_temporal.png
    """

    stub:
    """
    echo "gene,r_squared,slope,date_range" > ${gene_name}_temporal.csv
    echo "${gene_name},0.1,0.001,2015-2024" >> ${gene_name}_temporal.csv
    touch ${gene_name}_temporal.png
    """
}

/*
 * =============================================================================
 * Module 13: Variant Naming (CARD)
 * =============================================================================
 */

process VARIANT_NAMING {
    tag "${gene_name}"
    label 'process_low'
    publishDir "${params.outdir}/results/variants", mode: 'copy'

    input:
    tuple val(gene_name), path(fasta)

    output:
    tuple val(gene_name), path("${gene_name}_variants.csv"), emit: variants

    script:
    """
    card_lookup.py \
        --input ${fasta} \
        --gene ${gene_name} \
        --output ${gene_name}_variants.csv
    """

    stub:
    """
    echo "accession,gene,variant,identity" > ${gene_name}_variants.csv
    echo "STUB001,${gene_name},Unknown,100.0" >> ${gene_name}_variants.csv
    """
}

/*
 * =============================================================================
 * Module 18: Reporting
 * =============================================================================
 */

process REPORTING {
    tag "${gene_name}"
    label 'process_low'
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    tuple val(gene_name), path(qc_stats), path(aln_stats), path(saturation),
          path(selection_corrected), path(hgt_report), path(domains),
          path(temporal)

    output:
    tuple val(gene_name), path("${gene_name}_report.md"),   emit: report
    tuple val(gene_name), path("${gene_name}_report.html"), emit: report_html, optional: true

    script:
    """
    generate_report.py \
        --gene ${gene_name} \
        --qc-stats ${qc_stats} \
        --aln-stats ${aln_stats} \
        --saturation ${saturation} \
        --selection ${selection_corrected} \
        --hgt ${hgt_report} \
        --domains ${domains} \
        --temporal ${temporal} \
        --template ${projectDir}/templates/report.jinja2 \
        --out-md ${gene_name}_report.md \
        --out-html ${gene_name}_report.html
    """

    stub:
    """
    echo "# ${gene_name} Report (stub)" > ${gene_name}_report.md
    echo "<html>stub</html>" > ${gene_name}_report.html
    """
}
