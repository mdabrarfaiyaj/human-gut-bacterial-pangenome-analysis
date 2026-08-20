process BUILD_PHYLOGENY {
    tag "iqtree"
    publishDir "${params.outdir}/phylogeny", mode: 'copy'
    cpus 1
    memory '1.5 GB'

    input:
    path core_aln

    output:
    path "core_gene_alignment.aln.treefile", emit: tree
    path "core_gene_alignment.aln.*", emit: all_outputs

    script:
    def extra_args = task.ext.args ?: ''
    """
    iqtree2 \\
        -s ${core_aln} \\
        -m MFP \\
        -T ${task.cpus} \\
        ${extra_args}
    """
}
