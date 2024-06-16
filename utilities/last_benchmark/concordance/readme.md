### get data set
```
# https://github.com/genome-in-a-bottle/giab_data_indexes
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/HG004_HiSeq300x_fastq/140818_D00360_0047_BHA66FADXX/Project_RM8392/Sample_4A1/4A1_CGATGT_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/HG004_HiSeq300x_fastq/140818_D00360_0047_BHA66FADXX/Project_RM8392/Sample_4A1/4A1_CGATGT_L001_R2_001.fastq.gz

wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/NIST_Illumina_2x250bps/reads/D3_S1_L001_R1_001.fastq.gz	
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/NIST_Illumina_2x250bps/reads/D3_S1_L001_R2_001.fastq.gz

```


### Align BIO150
```
nthreads=12

workspace=/home/yiqing
r1=$workspace/yan/input/BIO150/4A1_CGATGT_L001_R1_001.fastq.gz
r2=$workspace/yan/input/BIO150/4A1_CGATGT_L001_R2_001.fastq.gz
output=$workspace/yan/output-bio150

#r1=$workspace/yan/input/BIO250/D3_S1_L001_R1_001.fastq.gz
#r2=$workspace/yan/input/BIO250/D3_S1_L001_R2_001.fastq.gz
#output=$workspace/yan/output-bio250


aligner=bwa
time $workspace/code/bwa/bwa mem -M -t $nthreads $workspace/data/bwa/whole_sequence.fna $r1 $r2 > $output/$aligner.sam

aligner=bowtie2
time $workspace/code/bowtie2/bowtie2 -p $nthreads -x $workspace/data/bowtie2/bowtie2 \
-1 $r1 -2 $r2 > $output/$aligner.sam

aligner=minimap2
time $workspace/code/minimap2/minimap2 -t $nthreads -ax sr $workspace/data/fsva-hg37/hg37.mmi \
 $r1 $r2 > $output/$aligner.sam

aligner=strobealign
time $workspace/code/strobealign/build/strobealign --use-index $workspace/data/fsva-hg37/hg37.fna \
-t $nthreads $r1 $r2  > $output/$aligner.sam

aligner=accalign
time $workspace/code/accel-align-release/accalign -I H -l 32 -t $nthreads -o $output/${aligner}.sam \
$workspace/data/fsva-hg37/hg37.fna $r1 $r2 

1>>$workspace/yan/log 2>&1

```