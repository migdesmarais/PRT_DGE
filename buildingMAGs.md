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


## Map reads to MAGs and ID who is present in diff treatments
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
 ```

## ID Taxonomy of unmapped reads

```
samtools addreplacerg -r "ID:S1\tSM:sample1" -o 55.rg.bam 55.sorted.bam
samtools addreplacerg -r "ID:S1\tSM:sample1" -o atm.rg.bam atm.sorted.bam

samtools merge merged_sorted_ATM_55.bam atm.rg.bam 55.rg.bam
samtools index merged_sorted_ATM_55.bam

samtools view -b -f 4 merged_sorted_ATM_55.bam > unmapped_ATM_55.bam
samtools fastq unmapped_ATM_55.bam -1 unmapped_R1.fastq -2 unmapped_R2.fastq -s unmapped_unpaired.fastq

cat unmapped_R1.fastq unmapped_R2.fastq unmapped_unpaired.fastq > all_unmapped.fastq

conda create -n kraken_env -c bioconda -c conda-forge kraken2
conda activate kraken_env
conda install -c bioconda krona

kraken2 --db /data_store/kraken_database \
  --threads 12 \
  --output kraken_output/kraken_all_unmapped.out \
  --report kraken_output/kraken_all_unmapped.report \
  --use-names \
  --single all_unmapped.fastq

bash ~/.conda/envs/kraken2_env/opt/krona/updateTaxonomy.sh


cut -f2,3 kraken_output/kraken_all_unmapped.out > kraken_output/for_krona.txt
ktImportTaxonomy -o kraken_output/unmapped_krona.html kraken_output/for_krona.txt

scp -r mdesmarais@fram.ucsd.edu:/scratch/mdesmarais/PRT_DGE/co-assembly/megahit_out/kraken_output/unmapped_krona.html ~/Downloads/
```
  
## Predict ORFs with prodigal
```
prodigal -i flavo.fasta -a prodigal/flavo_proteins.faa -d prodigal/flavo_nucleotides.fna -f gff -o prodigal/flavo_genes.gff -p meta
seqkit stats flavo_nucleotides.fna
seqkit stats flavo_proteins.faa

prodigal -i marini.fasta -a prodigal/marini_proteins.faa -d prodigal/marini_nucleotides.fna -f gff -o prodigal/marini_genes.gff -p meta
seqkit stats marini_nucleotides.fna
seqkit stats marini_proteins.faa

prodigal -i mori.fasta -a prodigal/mori_proteins.faa -d prodigal/mori_nucleotides.fna -f gff -o prodigal/mori_genes.gff -p meta
seqkit stats mori_nucleotides.fna
seqkit stats mori_proteins.faa

(seqkit) [mdesmarais@fram prodigal]$ seqkit stats marini_nucleotides.fna
file                    format  type  num_seqs    sum_len  min_len  avg_len  max_len
marini_nucleotides.fna  FASTA   DNA      3,708  3,996,234       69  1,077.7   12,237
(seqkit) [mdesmarais@fram prodigal]$ seqkit stats flavo_nucleotides.fna
file                   format  type  num_seqs    sum_len  min_len  avg_len  max_len
flavo_nucleotides.fna  FASTA   DNA      3,326  3,182,712       60    956.9    9,678
(seqkit) [mdesmarais@fram prodigal]$ seqkit stats flavo_proteins.faa
file                format  type     num_seqs    sum_len  min_len  avg_len  max_len
flavo_proteins.faa  FASTA   Protein     3,326  1,060,904       20      319    3,226
(seqkit) [mdesmarais@fram prodigal]$ seqkit stats marini_nucleotides.fna
file                    format  type  num_seqs    sum_len  min_len  avg_len  max_len
marini_nucleotides.fna  FASTA   DNA      3,708  3,996,234       69  1,077.7   12,237
(seqkit) [mdesmarais@fram prodigal]$ seqkit stats marini_proteins.faa
file                 format  type     num_seqs    sum_len  min_len  avg_len  max_len
marini_proteins.faa  FASTA   Protein     3,708  1,332,078       23    359.2    4,079
(seqkit) [mdesmarais@fram prodigal]$ seqkit stats mori_nucleotides.fna
file                  format  type  num_seqs    sum_len  min_len  avg_len  max_len
mori_nucleotides.fna  FASTA   DNA      3,978  3,924,498       66    986.6    8,484
(seqkit) [mdesmarais@fram prodigal]$ seqkit stats mori_proteins.faa
file               format  type     num_seqs    sum_len  min_len  avg_len  max_len
mori_proteins.faa  FASTA   Protein     3,978  1,308,166       22    328.9    2,828

```

## DRAM
```
git clone https://github.com/WrightonLabCSU/DRAM.git
cd DRAM
conda env create -f environment.yaml -n DRAM
conda activate DRAM
pip install .

seqkit seq -m 14 flavo_proteins.faa > flavo_proteins_14min.faa

DRAM.py annotate \
  -i flavo_proteins_min100.faa \
  -o DRAM/flavo_dram_output \
  --threads 8
```

## Eggnog mapper
```
conda create --name eggnog
conda activate eggnog
conda install -c bioconda -c conda-forge eggnog-mapper

mkdir /home/mdesmarais/.conda/envs/eggnog/lib/python2.7/site-packages/data
download_eggnog_data.py

emapper.py \
  -i flavo_proteins.faa \
  -o eggnog/flavo_emapper \
  --cpu 8 \
  -m diamond \
  --data_dir /home/mdesmarais/.conda/envs/eggnog/lib/python2.7/site-packages/data

emapper.py \
  -i mori_proteins.faa \
  -o eggnog/mori_emapper \
  --cpu 8 \
  -m diamond \
  --data_dir /home/mdesmarais/.conda/envs/eggnog/lib/python2.7/site-packages/data

emapper.py \
  -i marini_proteins.faa \
  -o eggnog/marini_emapper \
  --cpu 8 \
  -m diamond \
  --data_dir /home/mdesmarais/.conda/envs/eggnog/lib/python2.7/site-packages/data

```

## Kofamscan
```
conda create --name kofamscan

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda install kofamscan

mkdir kofam_db && cd kofam_db
wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
wget https://www.genome.jp/ftp/db/kofam/ko_list.gz
tar -xzf profiles.tar.gz
gunzip ko_list.gz



for hmm in /scratch/mdesmarais/PRT_DGE/MAGs/prodigal/kofam_db/profiles/*.hmm; do
  KO=$(basename "$hmm" .hmm)
  hmmsearch --cpu 1 --tblout kofamscan_manual/${KO}.tbl "$hmm" /scratch/mdesmarais/PRT_DGE/MAGs/prodigal/flavo_proteins.faa > /dev/null
done




exec_annotation \
  -o kofamscan/flavo_kofamscan.tsv \
  -f detail-tsv \
  -p /scratch/mdesmarais/PRT_DGE/MAGs/prodigal/kofam_db/profiles \
  -k /scratch/mdesmarais/PRT_DGE/MAGs/prodigal/kofam_db/ko_list \
  /scratch/mdesmarais/PRT_DGE/MAGs/prodigal/flavo_proteins.faa

exec_annotation \
  -o kofamscan/marini_kofamscan.tsv \
  -f detail-tsv \
  -p /scratch/mdesmarais/PRT_DGE/MAGs/prodigal/kofam_db/profiles \
  -k /scratch/mdesmarais/PRT_DGE/MAGs/prodigal/kofam_db/ko_list \
  /scratch/mdesmarais/PRT_DGE/MAGs/prodigal/marini_proteins.faa

exec_annotation \
  -o kofamscan/mori_kofamscan.tsv \
  -f detail-tsv \
  -p /scratch/mdesmarais/PRT_DGE/MAGs/prodigal/kofam_db/profiles \
  -k /scratch/mdesmarais/PRT_DGE/MAGs/prodigal/kofam_db/ko_list \
  /scratch/mdesmarais/PRT_DGE/MAGs/prodigal/mori_proteins.faa
```

## Prokka
```
conda create --name prokka
conda activate prokka
conda install -c conda-forge -c bioconda -c defaults prokka

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

interproscan.sh -i mori_proteins_clean.faa -f tsv -o output.tsv --goterms --iprlookup

```



