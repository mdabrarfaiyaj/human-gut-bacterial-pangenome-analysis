process BUILD_PANGENOME {
    tag "panaroo"
    publishDir "${params.outdir}/panaroo", mode: 'copy'
    cpus 1
    memory '1.5 GB'

    input:
    path gff_files

    output:
    path "gene_presence_absence.csv", emit: gene_presence_absence
    path "gene_presence_absence.Rtab", emit: rtab
    path "core_gene_alignment.aln", emit: core_aln
    path "panaroo_summary.txt", emit: summary

    script:
    """
    panaroo \\
        -i ${gff_files} \\
        -o . \\
        --clean-mode strict \\
        -a core \\
        --core_threshold ${params.core_threshold} \\
        -t ${task.cpus}

    # Panaroo writes a summary file; fall back to a generated one if absent
    if [ -f summary_statistics.txt ]; then
        cp summary_statistics.txt panaroo_summary.txt
    else
        {
          echo "Panaroo run summary"
          echo "core_threshold=${params.core_threshold}"
          echo "n_input_gff=\$(echo ${gff_files} | wc -w)"
        } > panaroo_summary.txt
    fi
    """
}
