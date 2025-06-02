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

