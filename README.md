# BLAST Pipeline and SNP Window Analysis

This document describes the workflow for blasting genes against a reference genome and finding genes near a SNP within a 1 Mbp window. All scripts and commands are designed for reproducibility and batch processing.

---

## Folder Structure

All files and scripts should be organized under:

```
Y:\Bigdata\computing\Shima\b_carinata\blast_gene_seed_color_gene_05_08_2025 (second_blast_with_TT_genes)
```

---

## Step 1: Create the BLAST Database

Generate a BLAST nucleotide database from your reference genome:

```sh
makeblastdb -in /home/mahmoudi/Bigdata/computing/Shima/b_carinata/blast_gene/GCA_040584065.1_ASM4058406v1_genomic.fna -dbtype nucl -out Brassica_car_DB
```

---

## Step 2: BLAST All Genes Against the Reference Genome

Blast each gene FASTA file against the reference database:

```sh
for gene in /home/mahmoudi/Bigdata/computing/Shima/b_carinata/blast_gene_seed_color_gene_05_08_2025/ar.gene/*.fa; do
  blastn \
    -query "$gene" \
    -db Brassica_car_DB \
    -out "${gene%.fa}_blast_results.txt" \
    -evalue 1e-5 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"
done
```

---

## Step 3: Merge All BLAST Results

Combine all BLAST results into a single file:

```sh
#!/bin/bash

# Output filename
output_file="all_genes_blast_results_with_qcov.tsv"

# Write the header
echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovs" > "$output_file"

# Append each *_blast_results.txt file
for file in /home/mahmoudi/Bigdata/computing/Shima/b_carinata/blast_gene_seed_color_gene_05_08_2025/ar.gene/*_blast_results.txt; do
  cat "$file" >> "$output_file"
done

echo "✅ All BLAST results (with query coverage) combined into: $output_file"
```

---

## Step 4: R Script to Find Genes Near a SNP (±1 Mbp)

Create a script named `find_genes_near_snp.R`:

```r
#!/usr/bin/env Rscript

# =========================
# Usage:
#   Rscript find_genes_near_snp.R CM081010.1 19070000
# Output: result_CM081010.1_19070000.txt
# =========================

suppressMessages({
  library(data.table)
})

# ---- Parse arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript find_genes_near_snp.R <chromosome> <position>")
}

chr <- args[1]
snp_pos <- as.numeric(args[2])

# ---- Parameters ----
window_size <- 1e6  # 1 Mbp

# ---- Output file name ----
out_file <- paste0("result_", chr, "_", snp_pos, ".txt")

# ---- Read BLAST results ----
blast_file <- "all_genes_blast_results_with_qcov.tsv"
blast <- fread(blast_file)

# Ensure columns are numeric
blast$sstart <- as.numeric(blast$sstart)
blast$send   <- as.numeric(blast$send)

# ---- Define window ----
min_pos <- snp_pos - window_size
max_pos <- snp_pos + window_size

# ---- Filter hits in ±1 Mbp ----
results <- blast[
  sseqid == chr &
  ((sstart >= min_pos & sstart <= max_pos) |
   (send   >= min_pos & send   <= max_pos))
]

# ---- Mark SNP_in_gene ----
if (nrow(results) > 0) {
  results[, SNP_in_gene := ifelse(
    (snp_pos >= pmin(sstart, send)) & (snp_pos <= pmax(sstart, send)),
    "YES", "NO"
  )]
  
  fwrite(results, file = out_file, sep = "\t")
} else {
  write("No gene copies found in ±1 Mbp window.", file = out_file)
}

message("Search completed. Results saved to: ", out_file)
```

---

## Notes

- Update all paths as necessary for your current environment.
- Ensure that BLAST+ and R with the `data.table` package are installed.
- This workflow will allow you to find any genes near a given SNP position within a ±1 Mbp window using the merged BLAST results.

---
