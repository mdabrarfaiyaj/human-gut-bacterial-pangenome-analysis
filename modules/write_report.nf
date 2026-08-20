process WRITE_REPORT {
    tag "report"
    publishDir "${params.outdir}/report", mode: 'copy'
    cpus 1
    memory '256 MB'

    input:
    path validation_report
    path qc_files
    path panaroo_summary
    path amr_files
    path figure_files

    output:
    path "final_report.md"

    script:
    """
    {
      echo "# Gut Bacterial Pangenome Project - Final Report"
      echo ""
      echo "Generated: \$(date)"
      echo ""
      echo "## Project question"
      echo ""
      echo "Which genes are conserved among strains of one gut bacterial species,"
      echo "which genes are variable, and are antimicrobial-resistance genes part"
      echo "of the accessory genome?"
      echo ""
      echo "## Validation"
      echo '```'
      cat ${validation_report}
      echo '```'
      echo ""
      echo "## Pangenome summary (Panaroo)"
      echo '```'
      cat ${panaroo_summary}
      echo '```'
      echo ""
      echo "## Figures"
      for f in ${figure_files}; do
        echo "- \$f"
      done
      echo ""
      echo "## Observation / Interpretation / Biological meaning / Limitation"
      echo ""
      echo "**Observation:** See figure_stats.txt and amr_heatmap.png for exact counts"
      echo "of core, accessory, and unique genes, and which strains carry which AMR genes."
      echo ""
      echo "**Interpretation:** Genes present in >=98% of genomes in this subset"
      echo "are treated as core; genes present in only some genomes are accessory."
      echo ""
      echo "**Biological meaning:** Variation in the accessory genome, including any"
      echo "AMR genes found there, may reflect strain-level functional differences"
      echo "acquired independently or via horizontal gene transfer."
      echo ""
      echo "**Limitation:** This analysis uses a 10-20 genome subset of the full"
      echo "Zou et al. 2019 Culturable Genome Reference and includes no phenotype"
      echo "data. Results describe genomic potential, not confirmed phenotypic"
      echo "resistance, and should not be generalized beyond this subset."
      echo ""
      echo "## Final claim"
      echo ""
      echo "> I performed a reproducible subset-level pangenome reanalysis of a"
      echo "> published gut bacterial genome resource (Zou et al. 2019)."
      echo ""
      echo "This is explicitly **not** a claim of reproducing the full 1,520-genome study."
    } > final_report.md
    """
}
