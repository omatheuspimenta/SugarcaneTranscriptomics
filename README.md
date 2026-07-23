# RNAseq to lncRNA in Sugarcane (Saccharum spp.)

This repository contains a step-by-step workflow to predict **long non-coding RNAs (lncRNAs)** from sugarcane RNA-seq data. This guide is written to be as reproducible and easy to understand as possible. 

---

## What are we doing here? (The Short Version)
1. **The Goal:** Find new, unmapped RNA sequences in sugarcane that do not code for proteins, are longer than 200 (base pairs), and might be long non-coding RNAs (lncRNAs).
2. **The Starting Point:** We start with data that has already been cleaned and aligned to a sugarcane reference genome (R570) using the `nf-core/rnaseq` pipeline.
3. **The Output:** A `FASTA` file (`lncRNA_candidates.fa`) containing the actual sequences of our candidate lncRNAs, ready for downstream analysis.

---

## Repository vs. Server Structure

Because genomic files are incredibly large (often hundreds of Gigabytes), only the code and lightweight tracking files are hosted here on GitHub. The actual raw data, BAM files, and large genome references remain on our servers.

**What is on GitHub:**
```text
.
├── README.md
├── code/
│   └── run_stringtie3.sh       # Script to run the assembly step
└── data/
    ├── lncRNA_candidate_ids.txt # Final list of candidate IDs
    ├── lncRNA_candidates.tar.xz # Compressed final FASTA results
    ├── references/              # Empty placeholder for reference genomes
    └── stringtie_denovo/
        └── mergelist_denovo.txt # List of all assembled transcripts

```

*Note: If you want to run this pipeline from scratch, you will need to download the R570 reference genome and place it in `data/reference/`, and place your aligned `.bam` files in the appropriate folders.*

---

## Prerequisites (Software Requirements)

To run this pipeline, you need a Unix/Linux environment with the following tools installed:

* **[StringTie](https://github.com/gpertea/stringtie)**: Assembles RNA reads into transcripts.
* **[gffcompare](https://github.com/gpertea/gffcompare)**: Compares our new transcripts against known genes.
* **[gffread](https://github.com/gpertea/gffread)**: Extracts DNA/RNA sequences based on genomic coordinates.

---

## Workflow

### 1. Set Up the References

Before starting, ensure your reference genome and annotations are in place.

* **Reference Annotation (GTF):** `data/reference/SofficinarumxspontaneumR570_771_v2.1.gene_exons.gtf`
* **Reference Genome (FASTA):** `data/reference/SofficinarumxspontaneumR570_771_v2_assembly.fasta`

### 2. Assemble Transcripts with StringTie

Using the aligned data (`.bam` files), we run `StringTie` to assemble the reads into continuous transcripts. We use a script (`run_stringtie3.sh`) to automate this for all samples.

The core command running inside the script looks like this:

```bash
stringtie "$bam_file" \
        -G data/reference/SofficinarumxspontaneumR570_771_v2.1.gene_exons.gtf \
        -o "$output_gtf" \
        -p 200 # Number of CPU threads

```

> **Note** 
> We run StringTie without the `-e` option to perform reference-guided transcriptome assembly, allowing the identification of novel transcripts and isoforms in addition to those already annotated. Using the `-e` option would restrict StringTie to quantifying only the transcripts present in the reference GTF, without assembling new transcripts.
 
After this runs for all samples, we create a file named `mergelist_denovo.txt` containing the paths to all the newly generated `.gtf` files.

### 3. Merge Transcripts into a Master Transcriptome

Since we have multiple samples (e.g., leaves, stems, different time points), we need to combine all their pieces into one single, unified master list of transcripts.

```bash
stringtie --merge \
  -G data/reference/SofficinarumxspontaneumR570_771_v2.1.gene_exons.gtf \
  -o data/sugarcane_master_transcriptome_denovo.gtf \
  data/stringtie_denovo/mergelist_denovo.txt

```

### 4. Isolate the "Unknown" Transcripts (gffcompare)

Now we have a giant list of transcripts. However, we only care about the new, unannotated ones. We use `gffcompare` to compare our master list against the known sugarcane genome. It will tag transcripts that overlap with known protein-coding genes so we can discard them.

```bash
gffcompare -r data/reference/SofficinarumxspontaneumR570_771_v2.1.gene_exons.gtf \
  -o gffcomp_denovo \
  data/sugarcane_master_transcriptome_denovo.gtf

```

> **Note**  
> The `-r` option specifies the reference annotation used for comparison. `gffcompare` assigns a class code to each assembled transcript, indicating whether it matches a known transcript exactly, represents a novel isoform, overlaps a known gene, or corresponds to a potentially novel intergenic transcript. This classification is commonly used to identify and filter novel transcript models.

### 5. Sanity Check

Let's see how many "new" transcripts we found. `gffcompare` assigns class codes:

* **`u`**: Unknown/intergenic (completely new, not overlapping any known gene).
* **`i`**: Intronic (falls entirely within the non-coding part of a known gene).
* **`x`**: Exonic overlap on the opposite strand.

Run this to count them:

```bash
cut -f3 gffcomp_denovo.sugarcane_master_transcriptome_denovo.gtf.tmap | sort | uniq -c | sort -rn

```

### 6. Filter by Size (> 200 base pairs)

By definition, **long** non-coding RNAs must be at least 200 nucleotides long. We use an `awk` command to filter our list, keeping only the IDs of transcripts that have tags `u`, `i`, or `x` AND are ≥ 200 bp.

```bash
awk -F'\t' 'NR>1 && ($3=="u" || $3=="x" || $3=="i") && $10>=200 { print $5 }' \
  gffcomp_denovo.sugarcane_master_transcriptome_denovo.gtf.tmap \
  > data/lncRNA_candidate_ids.txt

```

Check how many candidates we have:

```bash
wc -l data/lncRNA_candidate_ids.txt

```

### 7. Create a Filtered GTF File

A `GTF` file is just a list of coordinates. We now extract only the coordinates for our candidate lncRNAs from the master GTF.

```bash
awk -F'\t' 'NR==FNR { ids[$1]; next }
  {
    match($9, /transcript_id "([^"]+)"/, arr)
    if (arr[1] in ids) print
  }' data/lncRNA_candidate_ids.txt data/sugarcane_master_transcriptome_denovo.gtf \
  > data/lncRNA_candidates.gtf

```

### 8. Extract the Actual RNA Sequences (FASTA)

Downstream prediction softwares cannot read coordinates (`GTF`); they need the actual sequences in a `FASTA` format. We use `gffread` to fetch these letters from the reference genome based on our candidate coordinates.

```bash
gffread -w data/lncRNA_candidates.fa \
  -g data/reference/SofficinarumxspontaneumR570_771_v2_assembly.fasta \
  data/lncRNA_candidates.gtf

```

> **Note**
> The `-w` flag ensures we extract the **spliced** sequence, meaning the mature RNA sequence with the empty spaces/introns removed, which is exactly what prediction software expects).*

### 9. Final Validation

Let's make sure the number of sequences in our final FASTA file exactly matches the number of candidate IDs we filtered in step 6.

```bash
# Count sequences in the FASTA file
grep -c '^>' data/lncRNA_candidates.fa

# Count the original IDs
wc -l data/lncRNA_candidate_ids.txt

```

Both numbers should match!

