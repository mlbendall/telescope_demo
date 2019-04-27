# Telescope Demo


## 1. Setup conda environments

#####  Environment `singleloc_sim`

Contains the RNA-seq read simulator [polyester](https://www.bioconductor.org/packages/release/bioc/html/polyester.html)

```bash
conda create -n singleloc_sim bioconductor-polyester
```

#####  Environment `singleloc_align`

Alignment using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and default methods using [htseq-count](https://htseq.readthedocs.io/en/release_0.11.1/count.html)

```bash
conda create -n singleloc_align samtools bowtie2 htseq
```

#####  Environment `singleloc_tele`

Includes the [Telescope](https://github.com/mlbendall/telescope) package installed using pip.

```bash
conda create -n singleloc_tele \
  python=3.6 cython numpy=1.13 scipy=0.19.0 \
  intervaltree=2.1.0 pysam=0.12 samtools=1.5 openssl=1.0.2o

conda activate singleloc_tele
pip install git+git://github.com/mlbendall/telescope.git
conda deactivate
```

## 2. Setup references

Download human reference genome hg38 from ENCODE

```bash
wget -P refs https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
gunzip -c refs/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz > refs/hg38.fasta
```

Build index for bowtie2

```bash
conda activate singleloc_align
bowtie2-build refs/hg38.fasta refs/hg38
conda deactivate
```

Download HERV annotation

```bash
wget -O refs/HERV_rmsk.hg38.v2.2.gtf https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/HERV_rmsk.hg38.v2/transcripts.gtf 
```



## 3. Simulate reads

The script `simulate.R` reads FASTA sequences from [HML2_extracted.fna](refs/HML2_extracted.fna) and simulates using polyester.

```bash
conda activate singleloc_sim
Rscript simulate.R 1000 HML2_1q22
conda deactivate
```


## 4. Telescope approach

Align reads allowing multiple hits

```bash
bowtie2 \
 -k 100 --very-sensitive-local --score-min "L,0,1.6" \
 --rg-id sample_01 \
 -x refs/hg38 \
 -f -1 sims/sample_01_1.fasta -2 sims/sample_01_2.fasta \
 -S sims/multi.sam \
 2>&1 | tee sims/multi.log
```

```bash
conda activate singleloc_tele

telescope assign \
  --theta_prior 200000 \
  --max_iter 200 \
  --updated_sam \
  --outdir sims \
  sims/multi.sam \
  refs/HERV_rmsk.hg38.v2.gtf \
  2>&1 | tee sims/telescope.log

```

```bash
python tag_multi.py --reports refs/HERV_rmsk.hg38.v2.gtf < sims/multi.sam > sims/multi.tag.sam
cat sims/multi.tag.sam | grep -v 'ZT:Z:SEC' | samtools view -ub | samtools sort > sims/multi.tag.bam
samtools index sims/multi.tag.bam

perl -lane 'print if $F[3]==100;' < per_read.txt  | wc -l
perl -lane 'print if $F[2]>1;' < per_read.txt  | wc -l






# samtools view -ub -F 0x100 sims/telescope-updated.bam | samtools sort > sims/telescope-sorted.bam

samtools view -h sims/telescope-updated.bam | grep -v 'ZT:Z:SEC' | sed 's/YC:Z:209,236,228/YC:Z:255,255,255/' | samtools view -ub | samtools sort > sims/telescope-sorted.bam
samtools index sims/telescope-sorted.bam
```


# Default approach

```bash
bowtie2 \
 --rg-id sample_01 \
 -x refs/hg38 \
 -f -1 sims/sample_01_1.fasta -2 sims/sample_01_2.fasta \
 -S sims/default.sam \
 2>&1 | tee sims/default.log

samtools view -ub sims/default.sam | samtools sort > sims/default.bam
samtools index sims/default.bam

htseq-count -f bam -r pos -t exon -i locus sims/default.bam refs/HERV_rmsk.hg38.v2.gtf > sims/default.count1.tsv
htseq-count -f bam -r pos -t exon -i locus -a 0 sims/default.bam refs/HERV_rmsk.hg38.v2.gtf > sims/default.count0.tsv

grep -vP '\t0$' sims/default.count1.tsv
grep -vP '\t0$' sims/default.count0.tsv

# check the sum
# grep -vP '\t0$' sims/default.count0.tsv | cut -f2 | paste -sd+ | bc
```

Most reads are discarded due to low mapping quality.
If the mapping quality is ignored, expression is incorrectly detected from x HML2 loci.
Using our approach, we see that 100% of reads originating from HML2_1q22 locus are 
correctly assigned.


