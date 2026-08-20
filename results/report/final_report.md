# Gut Bacterial Pangenome Project - Final Report

Generated: Thu Aug 20 18:46:01 +06 2026

## Project question

Which genes are conserved among strains of one gut bacterial species,
which genes are variable, and are antimicrobial-resistance genes part
of the accessory genome?

## Validation
```
Validation report - Thu Aug 20 13:23:59 +06 2026
==============================

Genome files present : 16
Samples in metadata  : 16

Species check (must all be identical - one species per pangenome run):
Bacteroides vulgatus

PASS: single species, genome/metadata counts consistent.
```

## Pangenome summary (Panaroo)
```
Core genes	(99% <= strains <= 100%)	2572
Soft core genes	(95% <= strains < 99%)	0
Shell genes	(15% <= strains < 95%)	2816
Cloud genes	(0% <= strains < 15%)	3309
Total genes	(0% <= strains <= 100%)	8697```

## Figures
- amr_heatmap.png
- core_accessory_barplot.png
- phylogeny_amr.png

## Observation / Interpretation / Biological meaning / Limitation

**Observation:** See figure_stats.txt and amr_heatmap.png for exact counts
of core, accessory, and unique genes, and which strains carry which AMR genes.

**Interpretation:** Genes present in >=98% of genomes in this subset
are treated as core; genes present in only some genomes are accessory.

**Biological meaning:** Variation in the accessory genome, including any
AMR genes found there, may reflect strain-level functional differences
acquired independently or via horizontal gene transfer.

**Limitation:** This analysis uses a 16 genome subset of the full
Zou et al. 2019 Culturable Genome Reference and includes no phenotype
data. Results describe genomic potential, not confirmed phenotypic
resistance, and should not be generalized beyond this subset.

## Final claim

> I performed a reproducible subset-level pangenome reanalysis of a
> published gut bacterial genome resource (Zou et al. 2019).

This is explicitly **not** a claim of reproducing the full 1,520-genome study.
