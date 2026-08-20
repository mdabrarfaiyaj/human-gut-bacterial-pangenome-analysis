process VALIDATE_INPUT {
    tag "validate"
    publishDir "${params.outdir}/validation", mode: 'copy'

    input:
    path genomes
    path samplesheet

    output:
    path "validation_report.txt", emit: report

    script:
    """
    echo "Validation report - \$(date)" > validation_report.txt
    echo "==============================" >> validation_report.txt
    echo "" >> validation_report.txt

    n_genomes=\$(ls *.fna.gz 2>/dev/null | wc -l)
    n_samples=\$(tail -n +2 ${samplesheet} | wc -l)

    echo "Genome files present : \${n_genomes}" >> validation_report.txt
    echo "Samples in metadata  : \${n_samples}" >> validation_report.txt

    if [ "\${n_genomes}" -ne "\${n_samples}" ]; then
        echo "WARNING: genome file count does not match samplesheet row count" >> validation_report.txt
    fi

    echo "" >> validation_report.txt
    echo "Species check (must all be identical - one species per pangenome run):" >> validation_report.txt
    cut -f2 ${samplesheet} | tail -n +2 | sort -u >> validation_report.txt

    n_species=\$(cut -f2 ${samplesheet} | tail -n +2 | sort -u | wc -l)
    if [ "\${n_species}" -ne 1 ]; then
        echo "" >> validation_report.txt
        echo "ERROR: more than one species detected in samplesheet. A pangenome must not mix species." >> validation_report.txt
        cat validation_report.txt
        exit 1
    fi

    echo "" >> validation_report.txt
    echo "PASS: single species, genome/metadata counts consistent." >> validation_report.txt
    cat validation_report.txt
    """
}
