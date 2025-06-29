## Run pre-QC fastqc
```
conda activate fastqc
screen
fastqc *fastq.gz
```

##	Download and visualize QC files
```
scp -r mdesmarais@fram.ucsd.edu:/scratch/mdesmarais/PRT_DGE/raw_reads/NovaX25B_RNA/pre_fastqc/pre_fastqc ~/Downloads/
```

##	Calculate how many reads in raw files
```
echo -e "Filename\tReadCount" > read_counts.tsv

for f in *_R[12]_001.fastq.gz; do
    count=$(zcat "$f" | wc -l | awk '{print $1/4}')
    echo -e "$f\t$count" >> read_counts.tsv
done
```

##	Trimmomatic

#Make sample txt file
```
ls *_R1_001.fastq.gz | sed -E 's/_R1_001.fastq.gz//' > samples.txt
```
#run trimmomatic
```
conda activate trimmomatic_env

# Run Trimmomcatic
for file in $(cat samples.txt); do
  trimmomatic PE -phred33 -threads 12 \
    ${file}_R1_001.fastq.gz ${file}_R2_001.fastq.gz \
    /scratch/mdesmarais/PRT_DGE/reads_rna/${file}_paired_R1.fastq.gz \
    /scratch/mdesmarais/PRT_DGE/reads_rna/${file}_unpaired_R1.fastq.gz \
    /scratch/mdesmarais/PRT_DGE/reads_rna/${file}_paired_R2.fastq.gz \
    /scratch/mdesmarais/PRT_DGE/reads_rna/${file}_unpaired_R2.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
```

##	Calculate how many reads in trimmed read files
```
echo -e "Filename\tReadCount" > paired_read_counts.tsv

for f in *_paired_R[12].fastq.gz; do
    count=$(zcat "$f" | wc -l | awk '{print $1/4}')
    echo -e "$f\t$count" >> paired_read_counts.tsv
done
```

## Run post-QC fastqc
```
conda activate fastqc
screen
fastqc *fastq.gz
```

## Remove rRNA reads

```
conda create --name sortmerna_env
conda activate sortmerna_env
conda install sortmerna

# install database
wget https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz
mkdir rRNA_databases_v4
tar -xvf database.tar.gz -C rRNA_databases_v4

# run it
for R1 in *paired_R1.fastq.gz; do \
  R2=${R1/_R1/_R2}; \
  base=${R1%%_R1.fastq.gz}; \
  sortmerna \
    --ref /scratch/mdesmarais/PRT_DGE/reads_rna/sortmerna_db/rRNA_databases_v4/smr_v4.3_default_db.fasta \
    --reads $R1 --reads $R2 \
    --workdir rrna_sorted_reads/${base}_sortmerna \
    --aligned rrna_sorted_reads/${base}_rRNA \
    --other rrna_sorted_reads/${base}_non_rRNA \
    --fastx --log --threads 8; \
done
```

## Create transcriptome
# One way of doing things, assemble metatranscriptome and do a DGE in parallel with the MAGs way up approach.
```
conda create -n rna_spades_env python=3.9 spades=3.15.5 -c bioconda -c conda-forge
conda activate rna_spades_env

mkdir rnaSPAdes_output

for fq in *_paired_non_rRNA.fq.gz; do
    sample=$(basename "$fq" _paired_non_rRNA.fq.gz)
    echo "🔧 Assembling $sample..."

    rnaspades.py \
        -s "$fq" \
        -o "/scratch/mdesmarais/PRT_DGE/rnaSPAdes_output/$sample" \
        --threads 12 --memory 64
done


# ORFs with prodigal

```

## Map and count reads
```
conda create --name salmon_env salmon=1.10.1 -c bioconda -c conda-forge -y
conda activate salmon_env

salmon index -t flavo_nucleotides.fna -i salmon/flavo_nucleotides_index
salmon index -t marini_nucleotides.fna -i salmon/marini_nucleotides_index
salmon index -t mori_nucleotides.fna -i salmon/mori_nucleotides_index

nano run_salmon_quant.sh

#!/bin/bash

# Set directories
READ_DIR="/scratch/mdesmarais/PRT_DGE/rrna_sorted_reads"
INDEX_DIR="/scratch/mdesmarais/PRT_DGE/salmon_quant/salmon_indexes"
OUT_DIR="/scratch/mdesmarais/PRT_DGE/salmon_quant"

# Loop through each MAG index directory
for INDEX_PATH in "$INDEX_DIR"/*_index; do
    MAG_NAME=$(basename "$INDEX_PATH" _index)

    # Loop through each non-rRNA paired sample
    for R1 in "$READ_DIR"/*_paired_non_rRNA.fq.gz; do
        SAMPLE_NAME=$(basename "$R1" _paired_non_rRNA.fq.gz)
        SAMPLE_OUT="$OUT_DIR/${MAG_NAME}_${SAMPLE_NAME}"

        echo "Running Salmon for $SAMPLE_NAME on $MAG_NAME"

        mkdir -p "$SAMPLE_OUT"

        salmon quant -i "$INDEX_PATH" -l A \
            -r "$R1" \
            -o "$SAMPLE_OUT" \
            --validateMappings \
            --gcBias
    done
done

chmod +x run_salmon_quant.sh
conda activate salmon_env

# calculate hits
nano summarize_salmon_stats.sh

#!/bin/bash

# Directory containing Salmon quant outputs
SALMON_DIR="/scratch/mdesmarais/PRT_DGE/salmon_quant"

# Print header
echo -e "MAG\tSample\tTotal_Fragments\tMapped_Reads\tMapping_Rate(%)"

# Loop over each subdirectory
for D in "$SALMON_DIR"/*; do
    if [[ -d "$D" && -f "$D/lib_format_counts.json" ]]; then

        # Extract MAG and Sample info from folder name
        BASENAME=$(basename "$D")
        MAG=$(echo "$BASENAME" | cut -d'_' -f1)
        SAMPLE=$(echo "$BASENAME" | cut -d'_' -f2-)

        # Extract total fragments processed and mapped
        TOTAL=$(grep '"num_processed"' "$D/lib_format_counts.json" | grep -o '[0-9]\+')
        HITS=$(grep '"num_mapped"' "$D/lib_format_counts.json" | grep -o '[0-9]\+')

        # Compute mapping rate
        if [[ "$TOTAL" -gt 0 ]]; then
            RATE=$(awk "BEGIN {printf \"%.2f\", ($HITS/$TOTAL)*100}")
        else
            RATE="NA"
        fi

        # Output line
        echo -e "${MAG}\t${SAMPLE}\t${TOTAL}\t${HITS}\t${RATE}"
    fi
done

chmod +x summarize_salmon_stats.sh


```




