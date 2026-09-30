# RNA-seq Analysis Pipeline

This repository contains the bulk RNA-seq workflow used for post-alignment processing, quality control, gene-level quantification, and differential-expression analysis.

## Software versions

- STAR v2.5.2b
- SAMtools v1.16.1
- Subread / featureCounts v2.0.2
- Qualimap v2.3
- R v4.3.2
- DESeq2 v1.42.1
- data.table: version not recorded
- MultiQC: version not recorded

## Reference genome and annotation

Reads were aligned to the human GRCh38 reference genome.

Gene-level quantification used:

`Homo_sapiens.GRCh38.112.gtf`

The annotation path is provided as a command-line argument.

## Repository structure

```text
RNAseq_analysis_pipeline_v2/
├── README.md
├── software_versions.txt
├── scripts/
│   ├── 01_process_bams.sh
│   ├── 02_check_qc.sh
│   ├── 03_combine_counts.sh
│   ├── 04_multiqc_postalign.sh
│   └── 05_deseq2_covariate_adjusted.R
└── example/
    ├── bam_files.txt
    └── metadata_example.csv
```

## 1. BAM processing and gene-level quantification

Run:

```bash
bash scripts/01_process_bams.sh \
  --bam-list example/bam_files.txt \
  --gtf /path/to/Homo_sapiens.GRCh38.112.gtf \
  --qualimap /path/to/qualimap
```

For each BAM file, the script:

- creates a BAM index if needed;
- runs `samtools flagstat`;
- runs Qualimap BAM QC;
- runs featureCounts.

FeatureCounts parameters:

```text
-t exon
-g gene_id
-s 0
-p
-B
-C
```

Thus, gene counts are exon-based, summarized by `gene_id`, generated in unstranded paired-end mode, require both reads to be successfully aligned, and exclude chimeric fragments.

## 2. Sample-level QC

Run:

```bash
bash scripts/02_check_qc.sh example/bam_files.txt
```

Samples fail QC if either:

- mapped reads < 80%; or
- mean mapping quality < 30.

Outputs:

```text
qc_summary/passed_samples.txt
qc_summary/failed_samples.txt
qc_summary/qc_metrics.tsv
```

## 3. Combine featureCounts outputs

Create a list of count files:

```bash
find . -name "*_gene_counts.txt" | sort > gene_count_files.txt
```

Then run:

```bash
bash scripts/03_combine_counts.sh gene_count_files.txt combined_counts.csv
```

The script verifies that gene order is identical across samples and writes a gene-by-sample raw-count matrix.

## 4. MultiQC

After QC classification:

```bash
bash scripts/04_multiqc_postalign.sh
```

This generates MultiQC reports for all, passed, and failed samples when corresponding QC artifacts are available.

## 5. Covariate-adjusted DESeq2 analysis

`scripts/05_deseq2_covariate_adjusted.R` performs covariate-adjusted differential-expression analysis using DESeq2.

It uses the design:

```r
~ batch + AGE + Sex + Race + BigGroups
```

with `Control` as the reference group.

The script:

- removes genes with total counts <= 1;
- restricts analysis to samples present in both counts and metadata;
- excludes samples with missing design covariates;
- runs DESeq2;
- calculates:
  - Ischemic vs Control
  - Nonischemic vs Control
  - Nonischemic vs Ischemic
- writes normalized mean expression alongside DESeq2 results;
- generates a variance-stabilized matrix for QC/visualization.

Adjusted P values are those returned by DESeq2 using the Benjamini-Hochberg procedure.

## Input metadata format

The DESeq2 script expects:

```text
Sample,batch,AGE,Sex,Race,BigGroups
sample01,1,55,male,W,Control
sample02,1,60,female,B,Ischemic
```

Race is standardized to `W` or `B`. Sex is standardized to `male` or `female`.
