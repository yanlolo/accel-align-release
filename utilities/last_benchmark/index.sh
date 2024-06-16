# build AA index [based on master branch, elias branch not work]

workspace=/home/yiqing
ref=/home/yiqing/data/fsva-hg37/hg37.fna

## strobealign
time $workspace/code/strobealign/build/strobealign --create-index $ref -r 100

## minimap2
time $workspace/code/minimap2/minimap2 -d /home/yiqing/data/fsva-hg37/hg37.mmi $ref

## bwa
time $workspace/code/bwa/bwa index $ref

## bowtie2
time $workspace/code/bowtie2/bowtie2-build $ref /home/yiqing/data/fsva-hg37/hg37


## AA master branch
$workspace/code/accel-align-release/accindex -l 32 $workspace/data/fsva-hg37/hg37.fna
#mv hg37.fna.hash hg37.fna.kmer32.hash
#ln -s $workspace/data/fsva-hg37/hg37.fna.kmer32.hash ./hg37.fna.hash

## acc-strobmer index
time /media/ssd/ngs-data-analysis/code/accel-align-release/accindex \
--strobe-mode --create-index /media/ssd/ngs-data-analysis/data/fsva-hg37/hg37.fna -r 100




