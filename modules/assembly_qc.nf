process ASSEMBLY_QC {
    tag "${sample_id}"
    publishDir "${params.outdir}/qc", mode: 'copy'
    cpus 1
    memory '512 MB'

    input:
    tuple val(sample_id), path(fasta)

    output:
    path "${sample_id}.seqkit_stats.tsv", emit: stats

    script:
    """
    seqkit stats -a -T ${fasta} > ${sample_id}.seqkit_stats.tsv
    """
}
