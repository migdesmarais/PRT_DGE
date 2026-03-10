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
```

## Predict ORFs with prodigal
```
conda create --name prodigal
conda activate prodigal
conda install -c bioconda -c conda-forge prodigal

# Optional: activate environment
# conda activate prodigal_env

# Path to base directory
BASE_DIR="/scratch/mdesmarais/PRT_DGE/rnaSPAdes_output"
PRODIGAL_DIR="$BASE_DIR/prodigal_output"

# Create output folder if needed
mkdir -p "$PRODIGAL_DIR"

# Loop through all sample assemblies
for sample_dir in "$BASE_DIR"/*; do
    sample=$(basename "$sample_dir")
    
    # Skip the prodigal_output folder itself
    if [[ "$sample" == "prodigal_output" ]]; then
        continue
    fi

    fasta="$sample_dir/transcripts.fasta"
    outdir="$PRODIGAL_DIR/$sample"

    if [ -f "$fasta" ]; then
        echo "🧬 Running Prodigal for $sample"
        mkdir -p "$outdir"

        prodigal -i "$fasta" \
                 -a "$outdir/${sample}_proteins.faa" \
                 -d "$outdir/${sample}_nucleotides.fna" \
                 -f gff \
                 -o "$outdir/${sample}_genes.gff" \
                 -p meta

        # Optional summary
        if command -v seqkit &> /dev/null; then
            seqkit stats "$outdir/${sample}_nucleotides.fna"
            seqkit stats "$outdir/${sample}_proteins.faa"
        fi
    else
        echo "⚠️ transcripts.fasta not found for $sample"
    fi
done

seqkit stats *.fna
seqkit stats *.faa
```

## DRAM
```
# Ended up running DRAM on KBASE it's horrible to install
scp -r ~/Downloads/84_DRAM mdesmarais@fram.ucsd.edu:/scratch/mdesmarais/PRT_DGE/MAGs/prodigal/DRAM
scp -r ~/Downloads/117_DRAM mdesmarais@fram.ucsd.edu:/scratch/mdesmarais/PRT_DGE/MAGs/prodigal/DRAM
scp -r ~/Downloads/bin1_DRAM mdesmarais@fram.ucsd.edu:/scratch/mdesmarais/PRT_DGE/MAGs/prodigal/DRAM

git clone https://github.com/WrightonLabCSU/DRAM.git
cd DRAM
conda env create -f environment.yaml -n DRAM
conda activate DRAM
pip install .
cd /sc

DRAM-setup.py prepare_databases --output_dir DRAM_data --threads 8 &> dram_setup.log
DRAM-setup.py --output_dir DRAM_data --threads 8 &> dram_setup.log

seqkit seq -m 14 flavo_proteins.faa > flavo_proteins_14min.faa

DRAM.py annotate \
  -i flavo_proteins_min100.faa \
  -o DRAM/flavo_dram_output \
  --threads 
```

## Eggnog mapper
```
conda create --name eggnog
conda activate eggnog
conda install -c bioconda -c conda-forge eggnog-mapper

mkdir /home/mdesmarais/.conda/envs/eggnog/lib/python2.7/site-packages/data
cd /home/mdesmarais/.conda/envs/eggnog/lib/python2.7/site-packages/data
download_eggnog_data.py

# Define paths
INPUT_BASE="/scratch/mdesmarais/PRT_DGE/rnaSPAdes_output/prodigal_output"
OUTPUT_BASE="${INPUT_BASE}/eggnog"
DATA_DIR="/home/mdesmarais/.conda/envs/eggnog/lib/python2.7/site-packages/data"

# Make sure the eggnog output directory exists
mkdir -p "$OUTPUT_BASE"

# Loop through each sample directory
for sample_dir in "$INPUT_BASE"/*/; do
    sample=$(basename "$sample_dir")
    faa_file="${sample_dir}/${sample}_proteins.faa"
    output_prefix="${OUTPUT_BASE}/${sample}_emapper"

    if [ -f "$faa_file" ]; then
        echo "🔧 Annotating $sample..."

        emapper.py \
          -i "$faa_file" \
          -o "$output_prefix" \
          --cpu 8 \
          -m diamond \
          --data_dir "$DATA_DIR"
    else
        echo "⚠️ No protein file found for $sample, skipping."
    fi
done
```

## InterProScan
```
conda create --name interproscan
conda activate interproscan

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda install interproscan

mkdir my_interproscan
cd my_interproscan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.74-105.0/interproscan-5.74-105.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.74-105.0/interproscan-5.74-105.0-64-bit.tar.gz.md5

# Recommended checksum to confirm the download was successful:
md5sum -c interproscan-5.74-105.0-64-bit.tar.gz.md5
# Must return *interproscan-5.74-105.0-64-bit.tar.gz: OK*
# If not - try downloading the file again as it may be a corrupted copy.

tar -pxvzf interproscan-5.74-105.0-*-bit.tar.gz

python3 setup.py -f interproscan.properties

conda create -n java11-env openjdk=11 -y
conda activate java11-env
java -version

# Test with mobiD
./interproscan.sh -i test_all_appl.fasta -f tsv -dp -exclappl MobiDBLite
./interproscan.sh -i test_all_appl.fasta -f tsv

# Define base path (adjust if needed)
base="/scratch/mdesmarais/PRT_DGE/rnaSPAdes_output/prodigal_output"
output_dir="$base/interproscan"

# Make sure output dir exists
mkdir -p "$output_dir"

# Loop through each sample directory
for sample_dir in "$base"/*/; do
    sample=$(basename "$sample_dir")
    faa_file="${sample_dir}${sample}_proteins.faa"
    clean_faa="${sample_dir}${sample}_clean.faa"
    out_file="${output_dir}/${sample}_interpro.tsv"
    log_file="${output_dir}/${sample}_interpro.log"

    # Check if protein file exists
    if [[ -f "$faa_file" ]]; then
        echo "🔧 Cleaning $faa_file..."
        awk '/^>/ {match($0, /ID=([^;]+)/, a); print ">" a[1]; next} {gsub(/\*$/, "", $0); print}' \
            "$faa_file" > "$clean_faa"

        echo "🚀 Running InterProScan on $sample..."
        ./interproscan.sh \
            -i "$clean_faa" \
            -f tsv \
            -o "$out_file" \
            -exclappl MobiDBLite \
            -goterms &> "$log_file"
    else
        echo "⚠️  $faa_file not found, skipping $sample"
    fi
done

```
## Map reads to ORFs of transcriptome and count
```
conda create --name salmon_env salmon=1.10.1 -c bioconda -c conda-forge -y
conda activate salmon_env

# Loop through each sample directory and index fna files
base_dir="/scratch/mdesmarais/PRT_DGE/rnaSPAdes_output/prodigal_output"
index_dir="${base_dir}/salmon_index"
reads_dir="/scratch/mdesmarais/PRT_DGE/rrna_sorted_reads"
quant_dir="/scratch/mdesmarais/PRT_DGE/rnaSPAdes_output/prodigal_output/salmon_quant"

# Loop over all 6 samples for indexing
for sample in 551R_S399_L005 552R_S400_L005 553R_S401_L005 ATM1R_S402_L005 ATM2R_S403_L005 ATM3R_S404_L005; do
    fna_file="${base_dir}/${sample}/${sample}_nucleotides.fna"
    index_dir="${salmon_base}/${sample}"

    echo "📦 Indexing ${sample} into ${index_dir}"
    salmon index -t "$fna_file" -i "$index_dir" --threads 8
done

# QUANTIFYING
# List of sample names
samples=(551R_S399_L005 552R_S400_L005 553R_S401_L005 ATM1R_S402_L005 ATM2R_S403_L005 ATM3R_S404_L005)

# Loop over samples
for sample in "${samples[@]}"; do
  echo "🧬 Quantifying $sample..."

  salmon quant \
    -i "${index_base}/${sample}" \
    -l A \
    -r "${reads_dir}/${sample}_paired_non_rRNA.fq.gz" \
    -p 8 \
    --validateMappings \
    -o "${quant_dir}/${sample}_quant"
done

# we got about ~70% mapping rate

```

















## Map and count reads
## one way of doing thins is to annotate genome first then map reads to it, doing it here
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




