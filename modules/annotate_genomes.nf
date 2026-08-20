process ANNOTATE_GENOMES {
    tag "${sample_id}"
    publishDir "${params.outdir}/annotation/${sample_id}", mode: 'copy'
    cpus 1
    memory '1.5 GB'
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1

    input:
    tuple val(sample_id), path(fasta)
    output:
    tuple val(sample_id), path("${sample_id}.gff"), emit: gff
    path "${sample_id}*", emit: all_outputs
    script:
    """
    if [[ "${fasta}" == *.gz ]]; then
        gzip -d -c ${fasta} > uncompressed_genome.fasta
    else
        cp ${fasta} uncompressed_genome.fasta
    fi

    prokka \\
        --outdir prokka_out \\
        --prefix ${sample_id} \\
        --cpus ${task.cpus} \\
        --centre X --compliant \\
        uncompressed_genome.fasta
    cp prokka_out/${sample_id}.gff .
    cp prokka_out/${sample_id}.* . 2>/dev/null || true
    """
}
