/*
 * =============================================================================
 * Module 7: Gene Tree Building (IQ-TREE2)
 * =============================================================================
 * Maximum-likelihood phylogenetic tree with:
 *   - ModelFinder Plus (automatic model selection)
 *   - UFBoot (ultrafast bootstrap)
 *   - SH-aLRT support
 *   - Midpoint or outgroup rooting
 *
 * Nextflow concepts learned:
 *   - cpus and memory directives
 *   - label for resource allocation
 *   - Multiple output files from one tool
 * =============================================================================
 */

process BUILD_TREE {

    tag "${gene_name}"
    label 'process_high'

    publishDir "${params.outdir}/trees", mode: 'copy'

    input:
    tuple val(gene_name), path(alignment)

    output:
    tuple val(gene_name), path("${gene_name}.treefile"),  emit: tree
    tuple val(gene_name), path("${gene_name}.iqtree"),    emit: iqtree_log
    tuple val(gene_name), path("${gene_name}.log"),       emit: log

    script:
    """
    # Count sequences
    NSEQ=\$(grep -c "^>" ${alignment} || true)

    if [ "\$NSEQ" -ge 4 ]; then
        # Enough sequences for bootstrap
        iqtree \
            -s ${alignment} \
            --seqtype CODON11 \
            -m MFP \
            -bb 1000 \
            -alrt 1000 \
            -T AUTO \
            --threads-max ${task.cpus} \
            -pre ${gene_name}
    else
        # Too few for bootstrap, plain ML tree
        echo "[build_tree] Only \$NSEQ sequences, skipping bootstrap"
        iqtree \
            -s ${alignment} \
            --seqtype CODON11 \
            -m MFP \
            -T AUTO \
            --threads-max ${task.cpus} \
            -pre ${gene_name}
    fi

    # Root the tree (midpoint rooting)
    root_tree.py \
        --input ${gene_name}.treefile \
        --output ${gene_name}.treefile \
        --method midpoint
    """

    stub:
    """
    echo "((A:0.1,B:0.2):0.05,C:0.3);" > ${gene_name}.treefile
    echo "IQ-TREE stub log" > ${gene_name}.iqtree
    echo "stub log" > ${gene_name}.log
    """
}
