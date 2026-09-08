# miRNA-seq analysis in *Xenopus laevis*

Guide from raw FASTQ files to a count matrix and a list of
differentially expressed miRNAs using mirdeep2. 

Lines beginning with `$` are commands you type — don't type the `$` itself.
Lines beginning with `#` are comments. Everything else is output.

---

## Contents

1. [The dataset](#1-the-dataset)
2. [Install the software](#2-install-the-software)
3. [Set up the project folder](#3-set-up-the-project-folder)
4. [Download and format the miRNA reference](#4-download-and-format-the-mirna-reference)
5. [Download and format the genome](#5-download-and-format-the-genome)
6. [The raw data](#6-the-raw-data)
7. [Quality control](#7-quality-control)
8. [Map the reads with `mapper.pl`](#8-map-the-reads-with-mapperpl)
9. [Count the miRNAs with `quantifier.pl`](#9-count-the-mirnas-with-quantifierpl)
10. [Merge the samples into one matrix](#10-merge-the-samples-into-one-matrix)
11. [Write the metadata file](#11-write-the-metadata-file)
12. [Differential expression with DESeq2](#12-differential-expression-with-deseq2)
13. [Troubleshooting](#13-troubleshooting)

---

## 1. The dataset

NF stage 50 *X. laevis* tadpoles underwent spinal cord transection. Two tissues
were collected — neural progenitor/stem cells (NPSCs) and spinal cord fluid — at
2 days post-transection (dpt) and 6 dpt, plus uninjured controls.

|  | Uninjured | 2 dpt | 6 dpt |
|---|---|---|---|
| **NPSC** | 2 | 2 | 2 |
| **Fluid** | 2 | 2 | 2 |

Twelve libraries, `XS1` to `XS12`. Prepared with the QIAseq miRNA library kit and
sequenced on an Illumina NovaSeq 6000.

Two features of this design matter later:

- **Two replicates per group.** Enough to rank candidates, not enough to survive
  multiple-testing correction. Section 12 explains how we handle that.
- **An allotetraploid genome.** Some *X. laevis* miRNAs share identical mature
  sequences and cannot be told apart by short reads. They show identical counts.
  That's expected — see section 10.

---

## 2. Install the software

Ran using Unix/Linux. 
If on windows use Windows Subsystem for Linux (WSL): 
```bash
wsl -d Ubuntu
```


Install [Miniforge](https://github.com/conda-forge/miniforge) if you don't have
conda, then create one environment with everything in it:

```bash
$ conda create -n mirna -c conda-forge -c bioconda \
    fastqc multiqc mirdeep2 bowtie -y
```

Activate it. **Do this in every new terminal window:**

```bash
$ conda activate mirna
```

Check the software is installed — each should print a version or a help page:

```bash
$ fastqc --version
$ mapper.pl
$ bowtie --version
```

For the R side, install R (or RStudio), then run this once inside R:

```r
install.packages(c("BiocManager", "ggplot2"))
BiocManager::install("DESeq2")
```

---

## 3. Set up the project folder

```bash
$ mkdir -p ~/xla_mirna_analysis
$ cd ~/xla_mirna_analysis
```

We'll build this structure as we go:

```
~/xla_mirna_analysis/
├── raw_fastq_files/          your FASTQ files
├── mirna_ref/                miRBase mature + hairpin
├── genome_ref/               genome + bowtie index
├── qc/                       FastQC and MultiQC reports
└── mirdeep2/output/XS1..12/  one folder per sample
```

Every command below uses full paths starting with `~/`, so it doesn't matter
which folder you're in when you run it — except where a command writes into the
current folder, which is called out each time.

---

## 4. Download and format the miRNA reference

Two files: **mature** and **hairpin**.

```bash
$ mkdir -p ~/xla_mirna_analysis/mirna_ref
$ cd ~/xla_mirna_analysis/mirna_ref

$ wget https://mirbase.org/download/mature.fa
$ wget https://mirbase.org/download/hairpin.fa
```

### Keep only *X. laevis*

These files contain every species in miRBase. 

Take only X. laevis
```bash
$ awk '/^>/ {p = ($0 ~ /^>xla-/)} p' mature.fa  > mature_xla.fa
$ awk '/^>/ {p = ($0 ~ /^>xla-/)} p' hairpin.fa > hairpin_xla.fa
```

That awk reads as: at every header line, set a flag to whether it starts with
`>xla-`; print the line whenever the flag is on. 

### Remove whitespace from the IDs

miRBase headers carry a description after a space:

```
>xla-let-7a-5p MIMAT0046749 Xenopus laevis let-7a-5p
```

miRDeep2 refuses to run on IDs containing whitespace. It ships a script for this:

```bash
$ remove_white_space_in_id.pl mature_xla.fa  > mature_xla_renamed.fa
$ remove_white_space_in_id.pl hairpin_xla.fa > hairpin_xla_renamed.fa
```

Check both steps worked:

```bash
$ head -2 mature_xla_renamed.fa
$ grep -c ">" mature_xla_renamed.fa
```

You should see a bare `>xla-let-7a-5p` with nothing after it, and a few hundred
sequences.


---

## 5. Download and format the genome

The **X. laevis v10.1 assembly from Xenbase** (854 MB compressed):

```bash
$ mkdir -p ~/xla_mirna_analysis/genome_ref
$ cd ~/xla_mirna_analysis/genome_ref

$ wget https://download.xenbase.org/xenbase/Genomics/JGI/Xenla10.1/XENLA_10.1_genome.fa.gz
$ gunzip XENLA_10.1_genome.fa.gz
```

Same whitespace problem, same fix:

```bash
$ remove_white_space_in_id.pl XENLA_10.1_genome.fa > XENLA_10.1_genome_renamed.fa
```

### Build the bowtie index

miRDeep2 uses **bowtie 1**.

```bash
$ mkdir -p ~/xla_mirna_analysis/genome_ref/genome_index

$ bowtie-build XENLA_10.1_genome_renamed.fa \
    ~/xla_mirna_analysis/genome_ref/genome_index/bowtie_index_XENLA_10.1 \ --threads 16
```
change --threads based on computing power

You only do it once. 

It writes six files ending `.ebwt`. The thing you pass to `mapper.pl` in section 8
is the **prefix** — `.../bowtie_index_XENLA_10.1` — not any one of those files.

---

## 6. The raw data

```bash
$ cd ~/xla_mirna_analysis/raw_fastq_files
$ ls
```

```
XS10_S10_R1_001.fastq.gz  XS1_S1_R1_001.fastq.gz  XS4_S4_R1_001.fastq.gz  XS7_S7_R1_001.fastq.gz
XS11_S11_R1_001.fastq.gz  XS2_S2_R1_001.fastq.gz  XS5_S5_R1_001.fastq.gz  XS8_S8_R1_001.fastq.gz
XS12_S12_R1_001.fastq.gz  XS3_S3_R1_001.fastq.gz  XS6_S6_R1_001.fastq.gz  XS9_S9_R1_001.fastq.gz
```

Decompress if files gzipped, because **miRDeep2 cannot read gzipped files**:

```bash
$ gunzip -k *.gz
```

`-k` keeps the `.gz` originals. 

---

## 7. Quality control

```bash
$ mkdir -p ~/xla_mirna_analysis/qc

$ fastqc -t 4 -o ~/xla_mirna_analysis/qc \
    ~/xla_mirna_analysis/raw_fastq_files/*.fastq.gz

$ multiqc ~/xla_mirna_analysis/qc -o ~/xla_mirna_analysis/qc
```

Open `~/xla_mirna_analysis/qc/multiqc_report.html` in a browser.

---

## 8. Map the reads with `mapper.pl`

`mapper.pl` does four jobs at once: clips the adapter, discards reads shorter
than 18 nt, collapses identical reads into one entry with a count, and maps them
to the genome.

### The adapter

These are QIAseq miRNA libraries, so the 3′ adapter is `AACTGTAGGCACCATCAAT`.


### Run one sample

`mapper.pl` writes into whatever folder you're in, so make a folder per sample and
work inside it:

```bash
$ mkdir -p ~/xla_mirna_analysis/mirdeep2/output/XS1
$ cd ~/xla_mirna_analysis/mirdeep2/output/XS1

$ mapper.pl \
    ~/xla_mirna_analysis/raw_fastq_files/XS1_S1_R1_001.fastq \
    -e \
    -h \
    -j \
    -k AACTGTAGGCACCATCAAT \
    -l 18 \
    -m \
    -p ~/xla_mirna_analysis/genome_ref/genome_index/bowtie_index_XENLA_10.1 \
    -s reads_collapsed.fa \
    -t reads_collapsed_vs_genome.arf \
    -v \
    -o 4
```

```
# e = input is FASTQ
# h = convert to FASTA
# j = discard reads containing letters other than A, C, G, T, U, N
# k = 3' adapter sequence to clip
# l = discard reads shorter than 18 nt after clipping
# m = collapse identical reads to unique sequences with a count
# p = bowtie index prefix to map against
# s = write collapsed reads to this FASTA
# t = write mapping results to this ARF file
# v = print progress
# o = threads for bowtie
```


### Run the other eleven

Same command, one folder per sample.

```bash
$ for i in $(seq 1 12); do
    mkdir -p ~/xla_mirna_analysis/mirdeep2/output/XS${i}
    cd ~/xla_mirna_analysis/mirdeep2/output/XS${i}
    mapper.pl \
        ~/xla_mirna_analysis/raw_fastq_files/XS${i}_S${i}_R1_001.fastq \
        -e -h -j -k AACTGTAGGCACCATCAAT -l 18 -m \
        -p ~/xla_mirna_analysis/genome_ref/genome_index/bowtie_index_XENLA_10.1 \
        -s reads_collapsed.fa \
        -t reads_collapsed_vs_genome.arf \
        -v -o 4
done
```

---

## 9. Count the miRNAs with `quantifier.pl`

Still inside the sample's folder, so it picks up `reads_collapsed.fa` and writes
its output alongside:

```bash
$ cd ~/xla_mirna_analysis/mirdeep2/output/XS1

$ quantifier.pl \
    -p ~/xla_mirna_analysis/mirna_ref/hairpin_xla_renamed.fa \
    -m ~/xla_mirna_analysis/mirna_ref/mature_xla_renamed.fa \
    -r reads_collapsed.fa \
    -W \
    -d \
    -y XS1
```

```
# p = precursor (hairpin) sequences
# m = mature sequences
# r = collapsed reads from mapper.pl
# W = weight read counts by their number of mappings
# d = skip PDF generation (much faster)
# y = label appended to the output filenames
```


Output:

```
expression_XS1.html    miRNAs_expressed_all_samples_XS1.csv
expression_analyses    quantifier_run.log
```

The file we want is the `.csv`. Despite the extension it's **tab-separated**:

```bash
$ head miRNAs_expressed_all_samples_XS1.csv
```

```
#miRNA          read_count  precursor   total       seq         seq(norm)
xla-let-7a-5p   331691.23   xla-let-7a  331691.23   331691.23   24783.24
xla-let-7a-3p   612.00      xla-let-7a  612.00      612.00      45.73
xla-let-7b-5p   45633.42    xla-let-7b  45633.42    45633.42    3409.63
xla-let-7b-3p   5.00        xla-let-7b  5.00        5.00        0.37
```

- **Counts are fractional.** That's `-W`. A read mapping to three precursors
  contributes a third to each instead of a whole count to all three. 
- **The `seq(norm)` column is miRDeep2's own normalisation. Ignore it.** DESeq2
  normalises internally and needs raw counts.

### Run the other eleven

```bash
$ for i in $(seq 1 12); do
    cd ~/xla_mirna_analysis/mirdeep2/output/XS${i}
    quantifier.pl \
        -p ~/xla_mirna_analysis/mirna_ref/hairpin_xla_renamed.fa \
        -m ~/xla_mirna_analysis/mirna_ref/mature_xla_renamed.fa \
        -r reads_collapsed.fa \
        -W -d -y XS${i}
done
```

---

## 10. Merge the samples into one matrix

Each sample now has its own `.csv`. This script stitches them into one table. Run
it from the folder holding the sample folders:

```bash
$ cd ~/xla_mirna_analysis/mirdeep2/output
$ python3 merge_quantifier_outputs.py
```

```
XS1 -> ./XS1/miRNAs_expressed_all_samples_XS1.csv
XS2 -> ./XS2/miRNAs_expressed_all_samples_XS2.csv
...

Wrote mirna_count_matrix.tsv: 412 miRNAs x 12 samples
```

```bash
$ head -5 mirna_count_matrix.tsv
```

```
miRNA           XS1     XS2     XS3     XS4     XS5    XS6     XS7     XS8     XS9     XS10    XS11    XS12
xla-let-7a-3p   612     1204    490     665     25     1476    975     2391    1613    1808    2327    1116
xla-let-7a-5p   331691  414533  182107  284419  29207  629940  228788  183613  166384  149936  194241  316611
xla-let-7b-3p   5       34      14      20      0      33      69      203     89      139     89      36
xla-let-7b-5p   45633   67282   23823   41502   3412   93600   24208   43404   27164   36270   30679   39426
```

### What the script does

A miRNA made from several precursors gets **one row per precursor**, and `-W`
splits its reads between them. `xla-let-7a-5p` above happens to have one
precursor, but many don't. The script sums those rows so each miRNA appears once
with its full count.


> **Identical counts are fine.** Some *X. laevis* miRNAs share an identical mature
> sequence, so a read matching it is genuinely unassignable and they end up with
> the same counts in every sample. We leave them as separate rows. 

---

## 11. Write the metadata file

A plain table saying which sample is which. **Fill this in with your real sample
order — the values below are a placeholder.**

```bash
$ cd ~/xla_mirna_analysis
$ nano metadata.csv
```

```
sample,tissue,condition
XS1,npsc,uninjured
XS2,npsc,uninjured
XS3,npsc,dpt2
XS4,npsc,dpt2
XS5,npsc,dpt6
XS6,npsc,dpt6
XS7,fluid,uninjured
XS8,fluid,uninjured
XS9,fluid,dpt2
XS10,fluid,dpt2
XS11,fluid,dpt6
XS12,fluid,dpt6
```

The sample names must match the column headers of `mirna_count_matrix.tsv`
exactly, and the `tissue` values must match what section 12 asks for — `npsc` and
`fluid`, lower case. Don't put spaces after the commas.

> **Why `dpt2` and not `2dpt`.** R prefixes names that start with a digit, so
> `2dpt` silently becomes `X2dpt` in DESeq2's output and your contrast names stop
> matching what you typed.

---

## 12. Differential expression with DESeq2

Open R in the project folder, or open RStudio and run
`setwd("~/xla_mirna_analysis")`.

### Load the data

```r
library(DESeq2)

raw  <- read.delim("mirdeep2/output/mirna_count_matrix.tsv", check.names = FALSE)
meta <- read.csv("metadata.csv", row.names = 1)

meta$condition <- factor(meta$condition, levels = c("uninjured", "dpt2", "dpt6"))

counts_all <- as.matrix(raw[, rownames(meta)])
rownames(counts_all) <- raw$miRNA

stopifnot(identical(colnames(counts_all), rownames(meta)))
```

The order of `levels` matters: **the first is the reference**, so `uninjured` goes
first and fold changes are reported relative to it. The `stopifnot` catches a
column/metadata mismatch, which would otherwise make every result silently wrong.

### Fit one model per tissue

NPSCs and spinal cord fluid are different sample types with different library
characteristics. Fitting them together would force DESeq2 to estimate one
dispersion trend across both, which fits neither well.

```r
run_tissue <- function(tis) {
  m <- meta[meta$tissue == tis, ]
  m$condition <- droplevels(m$condition)
  cts <- counts_all[, rownames(m)]

  dds  <- DESeqDataSetFromMatrix(cts, m, design = ~ condition)
  keep <- rowSums(counts(dds) >= 10) >= 2      # >=10 counts in >=2 samples
  dds  <- dds[keep, ]
  message(tis, ": ", nrow(dds), " miRNAs kept")

  DESeq(dds)
}

dds_npsc  <- run_tissue("npsc")
dds_fluid <- run_tissue("fluid")
```

### Check the samples cluster sensibly

Do this **before** looking at any list of miRNAs.

```r
vsd <- varianceStabilizingTransformation(dds_npsc, blind = TRUE)
plotPCA(vsd, intgroup = "condition")
```

> Use `varianceStabilizingTransformation()`, not `vst()`. They're the same
> transformation, but `vst()` needs at least 1,000 features by default and will
> error on a few hundred miRNAs.

You want the two replicates of each timepoint sitting near each other, and the
timepoints separating. If a single replicate sits far from its partner, note it
before you read anything into the results — at n=2 there is no way to tell an
outlier from real variation statistically, so the plot is the check.

### Get the contrasts

All three comparisons come from the **same fitted model**. Don't re-run DESeq2 on
subsets — you'd lose the shared dispersion estimates, which is what makes a
two-replicate design analysable at all.

```r
res <- results(dds_npsc, contrast = c("condition", "dpt2", "uninjured"))
summary(res)

res_df <- as.data.frame(res)
res_df$miRNA <- rownames(res_df)
res_df <- res_df[order(res_df$pvalue), ]

sig <- subset(res_df, !is.na(pvalue) & pvalue < 0.05 & abs(log2FoldChange) > 1)
head(sig, 20)

write.csv(res_df, "npsc_dpt2_vs_uninjured_all.csv", row.names = FALSE)
write.csv(sig,    "npsc_dpt2_vs_uninjured_sig.csv", row.names = FALSE)
```

`contrast = c("condition", "dpt2", "uninjured")` reads as: in the variable
`condition`, compare `dpt2` against `uninjured`. A positive log2FoldChange means
higher at 2 dpt. Swap in `dpt6` for the other two comparisons, and use `dds_fluid`
for the other tissue — six comparisons in total.

`scripts/deseq2.R` runs all six and writes the volcano plots and heatmaps.

The filter uses **`pvalue`, not `padj`**. 

With two replicates per group, DESeq2 has almost no information about within-group
variability for any individual miRNA. The per-miRNA p-values are noisy. There is a high false-positive rate. This is exploratory. 


---

## 13. Troubleshooting

**`mapper.pl: command not found`**
The environment isn't active in this terminal. Run `conda activate mirna`.

**`Error: ... has not allowed whitespaces in its first identifier`**
A FASTA header still has a space. Re-run `remove_white_space_in_id.pl` on it.

**`quantifier.pl` gives an empty or all-zero matrix**
Check for a stray `-t xla` and remove it. Otherwise check the adapter.

**Almost nothing survives trimming; no 22 nt peak**
Wrong adapter. Confirm `AACTGTAGGCACCATCAAT` against the FastQC overrepresented
sequences table.

**`bowtie-build` errors, or `mapper.pl` can't find the index**
The `-p` argument is the index *prefix*, not a file. It must match exactly what
you gave as the second argument to `bowtie-build`.

**Counts are not whole numbers**
Expected — that's `-W`. The merge script rounds them.

**The same miRNA appears on several rows in the `.csv`**
Expected — one row per mature/precursor pair. The merge script sums them.

**Two miRNAs have identical counts in every sample**
Expected — they share an identical mature sequence. Leave them.

**`merge_quantifier_outputs.py` says "not found"**
Run it from the folder containing the `XS1`, `XS2`, … folders, or set `ROOT` at
the top of the script.

**`vst()` errors with "less than 'nsub' rows"**
Expected at this feature count. Use `varianceStabilizingTransformation()`.

---

## Further reading

- [miRBase](https://www.mirbase.org/)
- [Xenbase](https://www.xenbase.org/) — X. laevis v10.1 assembly
- [miRDeep2](https://github.com/rajewsky-lab/mirdeep2) — Friedländer *et al.*, 2012
- [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) — Love, Huber & Anders, 2014
