/*
 * =============================================================================
 * Module 9: Selection Analysis (HyPhy)
 * =============================================================================
 * Runs 7 HyPhy selection methods per gene (or per GARD partition).
 * Uses bacterial genetic code (NCBI Code 11).
 *
 * All methods check for minimum 4 sequences before running.
 * With <4 sequences, placeholder JSON is output.
 *
 * CRITICAL: Uses "Bacterial and Plant Plastid" genetic code, NOT "Universal".
 * CRITICAL: Runs AFTER GARD (Module 8) on partitioned alignments.
 * =============================================================================
 */

process SELECTION_SLAC {
    tag "${gene_name}"
    label 'process_medium'
    publishDir "${params.outdir}/results/selection", mode: 'copy'

    input:
    tuple val(gene_name), path(alignment), path(tree)

    output:
    tuple val(gene_name), val("SLAC"), path("${gene_name}_SLAC.json"), emit: result

    script:
    """
    NSEQ=\$(grep -c "^>" ${alignment} || true)
    if [ "\$NSEQ" -ge 4 ]; then
        hyphy slac \
            --alignment ${alignment} \
            --tree ${tree} \
            --code 11 \
            --output ${gene_name}_SLAC.json
    else
        echo "[SLAC] Skipping: only \$NSEQ sequences (need >=4)"
        echo '{"MLE": {"headers": [], "content": []}, "skipped": true}' > ${gene_name}_SLAC.json
    fi
    """

    stub:
    """
    echo '{"MLE": {"headers": [], "content": []}}' > ${gene_name}_SLAC.json
    """
}

process SELECTION_FEL {
    tag "${gene_name}"
    label 'process_medium'
    publishDir "${params.outdir}/results/selection", mode: 'copy'

    input:
    tuple val(gene_name), path(alignment), path(tree)

    output:
    tuple val(gene_name), val("FEL"), path("${gene_name}_FEL.json"), emit: result

    script:
    """
    NSEQ=\$(grep -c "^>" ${alignment} || true)
    if [ "\$NSEQ" -ge 4 ]; then
        hyphy fel \
            --alignment ${alignment} \
            --tree ${tree} \
            --code 11 \
            --output ${gene_name}_FEL.json
    else
        echo "[FEL] Skipping: only \$NSEQ sequences (need >=4)"
        echo '{"MLE": {"headers": [], "content": []}, "skipped": true}' > ${gene_name}_FEL.json
    fi
    """

    stub:
    """
    echo '{"MLE": {"headers": [], "content": []}}' > ${gene_name}_FEL.json
    """
}

process SELECTION_MEME {
    tag "${gene_name}"
    label 'process_medium'
    publishDir "${params.outdir}/results/selection", mode: 'copy'

    input:
    tuple val(gene_name), path(alignment), path(tree)

    output:
    tuple val(gene_name), val("MEME"), path("${gene_name}_MEME.json"), emit: result

    script:
    """
    NSEQ=\$(grep -c "^>" ${alignment} || true)
    if [ "\$NSEQ" -ge 4 ]; then
        hyphy meme \
            --alignment ${alignment} \
            --tree ${tree} \
            --code 11 \
            --output ${gene_name}_MEME.json
    else
        echo "[MEME] Skipping: only \$NSEQ sequences (need >=4)"
        echo '{"MLE": {"headers": [], "content": []}, "skipped": true}' > ${gene_name}_MEME.json
    fi
    """

    stub:
    """
    echo '{"MLE": {"headers": [], "content": []}}' > ${gene_name}_MEME.json
    """
}

process SELECTION_FUBAR {
    tag "${gene_name}"
    label 'process_medium'
    publishDir "${params.outdir}/results/selection", mode: 'copy'

    input:
    tuple val(gene_name), path(alignment), path(tree)

    output:
    tuple val(gene_name), val("FUBAR"), path("${gene_name}_FUBAR.json"), emit: result

    script:
    """
    NSEQ=\$(grep -c "^>" ${alignment} || true)
    if [ "\$NSEQ" -ge 4 ]; then
        hyphy fubar \
            --alignment ${alignment} \
            --tree ${tree} \
            --code 11 \
            --output ${gene_name}_FUBAR.json
    else
        echo "[FUBAR] Skipping: only \$NSEQ sequences (need >=4)"
        echo '{"MLE": {"headers": [], "content": []}, "skipped": true}' > ${gene_name}_FUBAR.json
    fi
    """

    stub:
    """
    echo '{"MLE": {"headers": [], "content": []}}' > ${gene_name}_FUBAR.json
    """
}

process SELECTION_BUSTED {
    tag "${gene_name}"
    label 'process_high'
    publishDir "${params.outdir}/results/selection", mode: 'copy'

    input:
    tuple val(gene_name), path(alignment), path(tree)

    output:
    tuple val(gene_name), val("BUSTED"), path("${gene_name}_BUSTED.json"), emit: result

    script:
    """
    NSEQ=\$(grep -c "^>" ${alignment} || true)
    if [ "\$NSEQ" -ge 4 ]; then
        hyphy busted \
            --alignment ${alignment} \
            --tree ${tree} \
            --code 11 \
            --output ${gene_name}_BUSTED.json
    else
        echo "[BUSTED] Skipping: only \$NSEQ sequences (need >=4)"
        echo '{"test results": {"p-value": 1.0}, "skipped": true}' > ${gene_name}_BUSTED.json
    fi
    """

    stub:
    """
    echo '{"test results": {"p-value": 0.5}}' > ${gene_name}_BUSTED.json
    """
}

process SELECTION_ABSREL {
    tag "${gene_name}"
    label 'process_high'
    publishDir "${params.outdir}/results/selection", mode: 'copy'

    input:
    tuple val(gene_name), path(alignment), path(tree)

    output:
    tuple val(gene_name), val("aBSREL"), path("${gene_name}_aBSREL.json"), emit: result

    script:
    """
    NSEQ=\$(grep -c "^>" ${alignment} || true)
    if [ "\$NSEQ" -ge 4 ]; then
        hyphy absrel \
            --alignment ${alignment} \
            --tree ${tree} \
            --code 11 \
            --output ${gene_name}_aBSREL.json
    else
        echo "[aBSREL] Skipping: only \$NSEQ sequences (need >=4)"
        echo '{"branch attributes": {}, "skipped": true}' > ${gene_name}_aBSREL.json
    fi
    """

    stub:
    """
    echo '{"branch attributes": {}}' > ${gene_name}_aBSREL.json
    """
}

process SELECTION_RELAX {
    tag "${gene_name}"
    label 'process_high'
    publishDir "${params.outdir}/results/selection", mode: 'copy'

    input:
    tuple val(gene_name), path(alignment), path(tree)

    output:
    tuple val(gene_name), val("RELAX"), path("${gene_name}_RELAX.json"), emit: result

    script:
    """
    NSEQ=\$(grep -c "^>" ${alignment} || true)
    if [ "\$NSEQ" -ge 4 ]; then
        hyphy relax \
            --alignment ${alignment} \
            --tree ${tree} \
            --code 11 \
            --test Foreground \
            --output ${gene_name}_RELAX.json
    else
        echo "[RELAX] Skipping: only \$NSEQ sequences (need >=4)"
        echo '{"test results": {"relaxation or intensification parameter": 1.0, "p-value": 1.0}, "skipped": true}' > ${gene_name}_RELAX.json
    fi
    """

    stub:
    """
    echo '{"test results": {"relaxation or intensification parameter": 1.0, "p-value": 0.5}}' > ${gene_name}_RELAX.json
    """
}
