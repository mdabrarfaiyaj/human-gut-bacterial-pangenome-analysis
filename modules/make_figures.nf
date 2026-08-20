process MAKE_FIGURES {
    tag "figures"
    publishDir "${params.outdir}/figures", mode: 'copy', overwrite: true
    publishDir "${params.figures_dir}", mode: 'copy', overwrite: true
    cpus 1
    memory '1 GB'

    input:
    path gene_presence_absence
    path amr_files
    path tree
    path qc_files

    output:
    path "*.png", emit: figures
    path "figure_stats.txt", emit: stats

    script:
    """
    mkdir -p amr_dir
    for f in ${amr_files}; do cp \$f amr_dir/; done

    ${params.python_bin} ${projectDir}/scripts/make_figures.py \\
        --gene-presence-absence ${gene_presence_absence} \\
        --amr-dir amr_dir \\
        --tree ${tree} \\
        --outdir .

    if [ ! -f phylogeny_amr.png ]; then
        echo "ERROR: phylogeny_amr.png was not generated" >&2
        exit 1
    fi
    """
}
