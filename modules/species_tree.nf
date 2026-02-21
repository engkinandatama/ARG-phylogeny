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

    # Build species tree
    iqtree2 \
        -s concat_housekeeping.fasta \
        --seqtype CODON11 \
        -m MFP \
        -bb 1000 \
        -alrt 1000 \
        -T AUTO \
        --threads-max ${task.cpus} \
        -pre species_ref

    # Root
    root_tree.py \
        --input species_ref.treefile \
        --output species_ref.treefile \
        --method midpoint
    """

    stub:
    """
    echo "((A:0.1,B:0.2):0.05,C:0.3);" > species_ref.treefile
    echo "stub" > species_ref.iqtree
    """
}
