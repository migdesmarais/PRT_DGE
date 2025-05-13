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

## Make directory in reads dir
```
mkdir -p megahit_assemblies
```

## Make script
```
nano run_megahit_all.sh
```

## Script:
```
for r1 in *_paired_R1.fastq.gz; do
  base=$(echo "$r1" | sed 's/_paired_R1.fastq.gz//')
  r2="${base}_paired_R2.fastq.gz"
  outdir="megahit_assemblies/${base}"

  echo "Assembling $base..."
  megahit \
    -1 "$r1" \
    -2 "$r2" \
    -o "$outdir" \
    --min-contig-len 1000 \
    --presets meta-sensitive \
    --num-cpu-threads 20
done
```

## Make it executable
```
chmod +x run_megahit_all.sh
```
## Run it
```
./run_megahit_all.sh
```

## Check assembly quality with quast
```
conda create --name quast_env
conda activate quast_env
conda install -c bioconda quast
quast --version
```

```
quast -o QC *_contigs.fa
scp -r migdesmarais@fram.ucsd.edu:/scratch/mdesmarais/OB_BONCAT-FACS-SEQ/reads/megahit_assemblies/contigs/QC /Users/migueldesmarais/Downloads
```

