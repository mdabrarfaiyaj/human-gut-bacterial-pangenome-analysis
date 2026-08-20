#!/usr/bin/env nextflow
/*
 * gut-pangenome: Gut Bacterial Pangenome Project
 *
 * Reproduces, at subset scale, the core-vs-accessory pangenome logic of
 * Zou et al. 2019 (Culturable Genome Reference) for ONE species
 * (Bacteroides vulgatus by default) using 16 genomes.
 *
 * Question: Which genes are conserved among strains of one gut bacterial
 * species, which genes are variable, and are antimicrobial-resistance
 * genes part of the accessory genome?
 *
 * Stages: VALIDATE_INPUT -> ASSEMBLY_QC -> ANNOTATE_GENOMES ->
 *         BUILD_PANGENOME -> DETECT_AMR -> BUILD_PHYLOGENY ->
 *         MAKE_FIGURES -> WRITE_REPORT
 */

nextflow.enable.dsl = 2

params.samplesheet   = "${projectDir}/data/metadata/samplesheet.tsv"
params.genome_dir    = "${projectDir}/data/genomes"
params.outdir        = "${projectDir}/results"
params.figures_dir   = "${projectDir}/figures"
params.core_threshold = 0.98

include { VALIDATE_INPUT   } from './modules/validate_input.nf'
include { ASSEMBLY_QC      } from './modules/assembly_qc.nf'
include { ANNOTATE_GENOMES } from './modules/annotate_genomes.nf'
include { BUILD_PANGENOME  } from './modules/build_pangenome.nf'
include { DETECT_AMR       } from './modules/detect_amr.nf'
include { BUILD_PHYLOGENY  } from './modules/build_phylogeny.nf'
include { MAKE_FIGURES     } from './modules/make_figures.nf'
include { WRITE_REPORT     } from './modules/write_report.nf'

workflow {

    // Channel of [sample_id, fasta_path] built from the samplesheet
    genomes_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def fasta = file("${params.genome_dir}/${row.sample_id}.fna.gz")
            tuple(row.sample_id, fasta)
        }

    VALIDATE_INPUT(genomes_ch.collect { it[1] }, file(params.samplesheet))

    ASSEMBLY_QC(genomes_ch)

    ANNOTATE_GENOMES(genomes_ch)

    BUILD_PANGENOME(ANNOTATE_GENOMES.out.gff.map { sample_id, gff -> gff }.collect())

    DETECT_AMR(genomes_ch)

    BUILD_PHYLOGENY(BUILD_PANGENOME.out.core_aln)

    MAKE_FIGURES(
        BUILD_PANGENOME.out.gene_presence_absence,
        DETECT_AMR.out.amr_tsv.collect(),
        BUILD_PHYLOGENY.out.tree,
        ASSEMBLY_QC.out.stats.collect()
    )

    WRITE_REPORT(
        VALIDATE_INPUT.out.report,
        ASSEMBLY_QC.out.stats.collect(),
        BUILD_PANGENOME.out.summary,
        DETECT_AMR.out.amr_tsv.collect(),
        MAKE_FIGURES.out.figures.collect()
    )
}

workflow.onComplete {
    log.info """
    ============================================
    gut-pangenome pipeline finished
    Status : ${workflow.success ? 'OK' : 'FAILED'}
    Results: ${params.outdir}
    Report : ${params.outdir}/report/final_report.md
    ============================================
    """
}
