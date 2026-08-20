process DETECT_AMR {
    tag "${sample_id}"
    publishDir "${params.outdir}/amr", mode: 'copy'
    cpus 1
    memory '1 GB'

    input:
    tuple val(sample_id), path(fasta)

    output:
    path "${sample_id}.amr.tsv", emit: amr_tsv

    script:
    """
    gunzip -c ${fasta} > tmp.fna

    amrfinder \\
        -n tmp.fna \\
        -o ${sample_id}.amr.tsv \\
        --name ${sample_id} \\
        --plus

    if [ ! -s ${sample_id}.amr.tsv ]; then
        echo "ERROR: amrfinder produced no output for ${sample_id}" >&2
        exit 1
    fi
    """
}
