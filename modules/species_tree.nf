/*
 * =============================================================================
 * Module 11: Species Tree (Housekeeping Genes)
 * =============================================================================
 * Builds a species reference tree from housekeeping genes by:
 *   1. Fetching housekeeping genes from NCBI
 *   2. Aligning each separately
 *   3. Concatenating alignments (common taxa)
 *   4. Building tree with IQ-TREE2
 *
 * Nextflow concepts learned:
 *   - Subworkflows (reusing fetch, align, tree modules)
 *   - .collect() to gather all gene alignments
 * =============================================================================
 */

process FETCH_HOUSEKEEPING {
    tag "${gene_name}"
    label 'process_internet'

    input:
    tuple val(gene_name), val(query), val(max_records), val(email)

    output:
    tuple val(gene_name), path("${gene_name}.fasta"), emit: fasta

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
    echo ">stub|gene=${gene_name}" > ${gene_name}.fasta
    echo "ATGATGATGATGATGATGATG" >> ${gene_name}.fasta
    echo "acc,gene,product,protein_id,country,date,length_nt" > ${gene_name}_meta.csv
    """
}

process ALIGN_HOUSEKEEPING {
    tag "${gene_name}"
    label 'process_medium'

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

process CONCAT_AND_BUILD_SPECIES_TREE {
    tag "species_tree"
    label 'process_high'

    publishDir "${params.outdir}/trees", mode: 'copy'

    input:
    path alignments  // All housekeeping gene alignments collected

    output:
    path "species_ref.treefile", emit: tree
    path "species_ref.iqtree",  emit: log

    script:
    """
    # Concatenate alignments by common taxa
    concat_alignments.py \
        --input-dir . \
        --output concat_housekeeping.fasta \
        --pattern "*.codon.fasta"

    # Check if concatenation produced a valid file
    if [ -s concat_housekeeping.fasta ]; then
        INPUT_ALN=concat_housekeeping.fasta
    else
        # Fallback: use the largest single HK gene alignment
        INPUT_ALN=\$(ls -S *.codon.fasta 2>/dev/null | head -1)
    fi

    if [ -n "\$INPUT_ALN" ] && [ -s "\$INPUT_ALN" ]; then
        NSEQ=\$(grep -c "^>" "\$INPUT_ALN" || true)
        echo "[species_tree] Using \$INPUT_ALN (\$NSEQ sequences)"

        if [ "\$NSEQ" -ge 4 ]; then
            iqtree \
                -s "\$INPUT_ALN" \
                --seqtype CODON11 \
                -m MFP \
                -bb 1000 \
                -alrt 1000 \
                -T AUTO \
                --threads-max ${task.cpus} \
                -pre species_ref
        elif [ "\$NSEQ" -ge 3 ]; then
            echo "[species_tree] Only \$NSEQ seqs, skipping bootstrap"
            iqtree \
                -s "\$INPUT_ALN" \
                --seqtype CODON11 \
                -m MFP \
                -T AUTO \
                --threads-max ${task.cpus} \
                -pre species_ref
        else
            echo "[species_tree] Only \$NSEQ seqs, creating placeholder tree"
            echo "((taxon_A:0.01,taxon_B:0.01):0.005,taxon_C:0.01);" > species_ref.treefile
            echo "Placeholder - too few sequences" > species_ref.iqtree
        fi

        # Root if treefile was created by iqtree
        if [ -f species_ref.treefile ] && grep -q ":" species_ref.treefile; then
            root_tree.py \
                --input species_ref.treefile \
                --output species_ref.treefile \
                --method midpoint
        fi
    else
        echo "[species_tree] WARNING: No usable alignments. Creating placeholder tree."
        echo "((taxon_A:0.01,taxon_B:0.01):0.005,taxon_C:0.01);" > species_ref.treefile
        echo "Placeholder tree - no alignments found" > species_ref.iqtree
    fi
    """

    stub:
    """
    echo "((A:0.1,B:0.2):0.05,C:0.3);" > species_ref.treefile
    echo "stub" > species_ref.iqtree
    """
}
