## Run pre-QC fastqc
```
conda activate fastqc
screen
fastqc *fastq.gz
```

##	Download and visualize QC files
```
scp -r mdesmarais@fram.ucsd.edu:/scratch/mdesmarais/PRT_DGE/raw_reads/NovaX25B_DNA/pre_fastqc/pre_fastqc ~/Downloads/
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

```
conda activate trimmomatic_env

# Run Trimmomcatic
for file in $(cat samples.txt); do
  trimmomatic PE -phred33 -threads 12 \
    ${file}_R1_001.fastq.gz ${file}_R2_001.fastq.gz \
    /scratch/mdesmarais/PRT_DGE/reads_dna/${file}_paired_R1.fastq.gz \
    /scratch/mdesmarais/PRT_DGE/reads_dna/${file}_unpaired_R1.fastq.gz \
    /scratch/mdesmarais/PRT_DGE/reads_dna/${file}_paired_R2.fastq.gz \
    /scratch/mdesmarais/PRT_DGE/reads_dna/${file}_unpaired_R2.fastq.gz \
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

##	Download and visualize QC files
```
scp -r mdesmarais@fram.ucsd.edu:/scratch/mdesmarais/PRT_DGE/reads_dna/post_fastqc ~/Downloads/
```

## Install MEGAHIT
```
conda create --name megahit
conda activate megahit
conda install -c bioconda megahit
megahit --version
```

## Install metaspades
```
conda create -n spades_env python=3.8 -y
conda activate spades_env
conda install -c bioconda spades
```

## Concatenate forward and reverse clean reads for co-assembly:
> `cat *_R1.fastq.gz > ALL_READS_1.fastq.gz`
> `cat *_R2.fastq.gz > ALL_READS_2.fastq.gz`


## Run metaspades on merged reads
```
metaspades.py \
  -1 /scratch/mdesmarais/PRT_DGE/co-assembly/ALL_READS_1.fastq.gz \
  -2 /scratch/mdesmarais/PRT_DGE/co-assembly/ALL_READS_2.fastq.gz \
  -o /scratch/mdesmarais/PRT_DGE/co-assembly/metaspades_out \
  -t 12 -m 128
```

## Run megahit on merged reads
```
megahit \
  -1 ALL_READS_1.fastq.gz \
  -2 ALL_READS_2.fastq.gz \
  -o megahit_out \
  -t 12 \
  -m 0.95
```

## QC on assemblies
```
quast \
  metaspades_out/contigs.fasta \
  megahit_out/final.contigs.fa \
  -1 ALL_READS_1.fastq.gz \
  -2 ALL_READS_2.fastq.gz \
  -o quast_comparison \
  -t 12
```

## Install binning env
```
conda create -n binning_tools_env -c bioconda -c conda-forge \
  metabat2 maxbin2 concoct bowtie2 samtools bbmap -y
```

## Install refine MAGs env
```
conda create -n refine_mags_env -c bioconda -c conda-forge \
  das_tool drep checkm-genome diamond prodigal -y
```

## Map original reads to final assembly - MEGAHIT was a lot better so going with the MEGAHIT assembly. See QC results.
```
# Build index
bowtie2-build final.contigs.fa assembly_index

# Map reads (for each sample)
bowtie2 -x assembly_index -1 /scratch/mdesmarais/PRT_DGE/reads_dna/ATMD_S501_L006_paired_R1.fastq.gz -2 /scratch/mdesmarais/PRT_DGE/reads_dna/ATMD_S501_L006_paired_R2.fastq.gz | samtools view -bS - | samtools sort -o ATMD.bam

bowtie2 -x assembly_index -1 /scratch/mdesmarais/PRT_DGE/reads_dna/55D_S500_L006_paired_R1.fastq.gz -2 /scratch/mdesmarais/PRT_DGE/reads_dna/55D_S500_L006_paired_R2.fastq.gz | samtools view -bS - | samtools sort -o 55D.bam

# Index BAM
samtools index ATMD.bam
samtools index 55D.bam

# Repeat for all samples (e.g., sample2.bam, etc.)
```

## Calculate contig depth
```
jgi_summarize_bam_contig_depths --outputDepth depth.txt ATMD.bam 55D.bam
```

## Run binning tools
```
# MetaBAT2
metabat2 -i final.contigs.fa -a depth.txt -o bins_metabat2/bin

# MaxBin2
run_MaxBin.pl -contig final.contigs.fa -abund depth.txt -out bins_maxbin2/bin

# CONCOCT (requires special input format)
cut_up_fasta.py final.contigs.fa -c 10000 -o 0 --merge_last -b concoct_input/contigs_10k.fa
concoct_coverage_table.py concoct_input/contigs_10k.fa ATMD.bam 55D.bam > concoct_input/coverage_table.tsv

concoct --composition_file concoct_input/contigs_10k.fa \
        --coverage_file concoct_input/coverage_table.tsv \
        -b concoct_output/

cut_up_fasta.py final.contigs.fa -c 10000 -o 0 --merge_last -b concoct_input/contigs_10K.bed > concoct_input/contigs_10K.fa
concoct_coverage_table.py concoct_input/contigs_10K.bed ATMD.bam 55D.bam > coverage_table.tsv
concoct --composition_file concoct_input/contigs_10K.fa --coverage_file concoct_input/coverage_table.tsv -b concoct_output/
merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv
mkdir concoct_output/fasta_bins

extract_fasta_bins.py final.contigs.clean.fa concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins
```

## Derep with DAStool

```
for f in /scratch/mdesmarais/PRT_DGE/co-assembly/megahit_out/bins_metabat2/*.fa; do
  bin=$(basename "$f" .fa)
  grep "^>" "$f" | sed 's/^>//' | awk -v b="$bin" '{print $0 "\t" b}' 
done > metabat2_contig2bin.tsv

for f in /scratch/mdesmarais/PRT_DGE/co-assembly/megahit_out/bins_maxbin2/*.fa; do
  bin=$(basename "$f" .fa)
  grep "^>" "$f" | sed 's/^>//' | awk -v b="$bin" '{print $0 "\t" b}' 
done > maxbin2_contig2bin.tsv

for f in /scratch/mdesmarais/PRT_DGE/co-assembly/megahit_out/bins_concoct/*.fa; do
  bin=$(basename "$f" .fa)
  grep "^>" "$f" | sed 's/^>//' | awk -v b="$bin" '{print $0 "\t" b}' 
done > concoct_contig2bin.tsv
```

```
DAS_Tool \
  -i /scratch/mdesmarais/PRT_DGE/co-assembly/megahit_out/metabat2_contig2bin.tsv,\
/scratch/mdesmarais/PRT_DGE/co-assembly/megahit_out/maxbin2_contig2bin.tsv,\
/scratch/mdesmarais/PRT_DGE/co-assembly/megahit_out/concoct_contig2bin.tsv \
  -l metabat2,maxbin2,concoct \
  -c /scratch/mdesmarais/PRT_DGE/co-assembly/megahit_out/final.contigs.clean.fa \
  -o /scratch/mdesmarais/PRT_DGE/co-assembly/megahit_out/dastool_output \
  --search_engine diamond \
  --write_bins \
  --threads 12 \
  --write_bin_evals

scp mdesmarais@fram.ucsd.edu:/scratch/mdesmarais/PRT_DGE/co-assembly/megahit_out/dastool_output_DASTool_summary.tsv ~/Downloads/
scp mdesmarais@fram.ucsd.edu:/scratch/mdesmarais/PRT_DGE/co-assembly/megahit_out/dastool_output_DASTool__bins ~/Downloads/
```

## Assign taxonomy
```
conda create -n gtdbtk_env -c bioconda -c conda-forge gtdbtk
conda activate gtdbtk_env
gtdbtk --help

export GTDBTK_DATA_PATH=/data_store/gtdbtk_db/release226

gtdbtk classify_wf \
  --genome_dir /scratch/mdesmarais/PRT_DGE/co-assembly/megahit_out/dastool_output_DASTool_bins/gtdbtk_ready_MAGs \
  --out_dir gtdbk_output \
  --cpus 12 \
  --extension fa \
  --mash_db PRT_DGE_mash.msh

gtdbtk classify_wf \
  --genome_dir gtdbtk_ready_MAGs \
  --out_dir gtdbtk_output \
  --cpus 12 \
  --mash_db your_custom_or_ref.msh

```

```
conda create -n kraken2_env -c bioconda -c conda-forge python=3.10 kraken2
conda activate kraken2_env
```




cat dastool_output_DASTool_bins/*.fa > all_MAGs.fna

bowtie2-build all_MAGs.fna all_MAGs_index

bowtie2 -x all_MAGs_index \
  -1 /scratch/mdesmarais/PRT_DGE/reads_dna/55D_S500_L006_paired_R1.fastq.gz -2 /scratch/mdesmarais/PRT_DGE/reads_dna/55D_S500_L006_paired_R2.fastq.gz \
  -S 55.sam -p 12

bowtie2 -x all_MAGs_index \
  -1 /scratch/mdesmarais/PRT_DGE/reads_dna/ATMD_S501_L006_paired_R1.fastq.gz -2 /scratch/mdesmarais/PRT_DGE/reads_dna/ATMD_S501_L006_paired_R2.fastq.gz \
  -S atm.sam -p 12

samtools view -bS 55.sam | samtools sort -o 55.sorted.bam
samtools index 55.sorted.bam

samtools view -bS atm.sam | samtools sort -o atm.sorted.bam
samtools index atm.sorted.bam

conda create -n coverm_env -c bioconda -c conda-forge coverm
conda activate coverm_env

for f in *.fa; do mv "$f" "${f%.fa}.fna"; done

coverm genome \
  -b atm.sorted.bam 55.sorted.bam \
  --genome-fasta-directory /scratch/mdesmarais/PRT_DGE/co-assembly/megahit_out/dastool_output_DASTool_bins \
  -m relative_abundance tpm rpkm \
  --min-read-aligned-percent 75 \
  --min-read-percent-identity 95 \
  --threads 12 \
  > coverm_results.tsv
 

conda create -n kraken_env -c bioconda -c conda-forge kraken2
conda activate kraken_env

mkdir -p kraken_gtdb_db
kraken2-build --download-taxonomy --db kraken_gtdb_db

for genome in /scratch/mdesmarais/PRT_DGE/co-assembly/megahit_out/dastool_output_DASTool_bins/*.fna; do
  kraken2-build --add-to-library "$genome" --db kraken_gtdb_db
done

kraken2-build --add-to-library /scratch/mdesmarais/PRT_DGE/co-assembly/megahit_out/dastool_output_DASTool_bins/*.fna --db /kraken_gtdb_db

kraken2-build --build --db kraken_gtdb_db

kraken2 --db ~/kraken_gtdb_db --paired sample_R1.fastq.gz sample_R2.fastq.gz \
  --report kraken_report.txt --output kraken_output.txt






