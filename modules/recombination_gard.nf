/*
 * =============================================================================
 * Module 8: Recombination Detection (GARD)
 * =============================================================================
 * Detects recombination breakpoints using HyPhy GARD.
 * Splits alignment into partitions for downstream selection analysis.
 *
 * CRITICAL: This MUST complete BEFORE selection analysis (Module 9).
 * The channel dependency in main.nf enforces this ordering.
 *
 * Nextflow concepts learned:
 *   - Strict process ordering via channel dependencies
 *   - Conditional outputs
 * =============================================================================
 */

process RECOMBINATION_GARD {

    tag "${gene_name}"
    label 'process_high'

    publishDir "${params.outdir}/results/recombination", mode: 'copy'

    input:
    tuple val(gene_name), path(alignment)

    output:
    tuple val(gene_name), path("${gene_name}_GARD.json"),        emit: gard_json
    tuple val(gene_name), path("${gene_name}_partitions/*.fasta"), emit: partitions, optional: true
    tuple val(gene_name), path(alignment),                        emit: original_alignment

    script:
    """
    NSEQ=\$(grep -c "^>" ${alignment} || true)

    if [ "\$NSEQ" -ge 4 ]; then
        # Run GARD (may crash on low-diversity data)
        hyphy gard \
            --alignment ${alignment} \
            --srv Yes \
            --output ${gene_name}_GARD.json \
        || {
            echo "[gard] GARD crashed, creating fallback (no recombination assumed)"
            echo '{"breakpoints": [], "skipped": true, "reason": "gard_error"}' > ${gene_name}_GARD.json
        }

        # Parse GARD results and split alignment into partitions
        mkdir -p ${gene_name}_partitions
        if python3 -c "import json; d=json.load(open('${gene_name}_GARD.json')); assert not d.get('skipped')" 2>/dev/null; then
            partition_gard.py \
                --gard-json ${gene_name}_GARD.json \
                --alignment ${alignment} \
                --gene ${gene_name} \
                --outdir ${gene_name}_partitions
        else
            cp ${alignment} ${gene_name}_partitions/${gene_name}_partition_1.fasta
        fi
    else
        echo "[gard] Only \$NSEQ sequences, skipping GARD (need >=4)"
        echo '{"breakpoints": [], "skipped": true, "reason": "too_few_sequences"}' > ${gene_name}_GARD.json
        mkdir -p ${gene_name}_partitions
        cp ${alignment} ${gene_name}_partitions/${gene_name}_partition_1.fasta
    fi
    """

    stub:
    """
    echo '{"breakpoints": []}' > ${gene_name}_GARD.json
    mkdir -p ${gene_name}_partitions
    cp ${alignment} ${gene_name}_partitions/${gene_name}_partition_1.fasta
    """
}
