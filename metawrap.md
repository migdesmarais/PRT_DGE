## Try metaWRAP
```
conda create --name metawrap
conda activate metawrap
conda install -c bioconda metawrap
metawrap --version
```

#### Concatenate forward and reverse clean reads for co-assembly:
> `cat *_R1.fastq.gz > ALL_READS_1.fastq.gz`
> `cat *_R2.fastq.gz > ALL_READS_2.fastq.gz`

### Assembly
Assemble reads into contigs with --megahit or --metaspades and create assembly report with QUAST
> ` metawrap assembly -1 ALL_READS_1.fastq -2 ALL_READS/ALL_READS_2.fastq -m 1024 -t 12 --megahit -o ASSEMBLY`
