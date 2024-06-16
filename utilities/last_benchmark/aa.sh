code_dir=/home/yiqing/code/
in_dir=/home/yiqing/yan/input/simulate/1m/
out_dir=/home/yiqing/yan/output
ref=/home/yiqing/data/fsva-hg37/hg37.fna
nthreads=12
rlen=100
r=sv-1m-${rlen}-r.fastq
truth=sv-1m-${rlen}-se-align.sam

## SE ##
aligner=accalign
time $code_dir/accel-align-release/accalign -I H -l 32 -t $nthreads -o $out_dir/${aligner}.sam $ref $in_dir/$r
python $code_dir/accel-align-release/utilities/match.py $in_dir/$truth $out_dir/${aligner}.sam $out_dir/${aligner}.csv

## PE ##
r1=sv-1m-${rlen}-r1.fastq
r2=sv-1m-${rlen}-r2.fastq
truth=sv-1m-${rlen}-pe-align.sam
aligner=accalign
time $code_dir/accel-align-release/accalign -I H -l 32 -t $nthreads -o $out_dir/${aligner}.sam $ref $in_dir/$r1 $in_dir/$r2
python $code_dir/accel-align-release/utilities/match.py $in_dir/$truth $out_dir/${aligner}.sam $out_dir/${aligner}.csv


## concor ##
nthreads=12

workspace=/home/yiqing
r1=$workspace/yan/input/BIO150/4A1_CGATGT_L001_R1_001.fastq.gz
r2=$workspace/yan/input/BIO150/4A1_CGATGT_L001_R2_001.fastq.gz
output=$workspace/yan/output-bio150
aligner=accalign
time $workspace/code/accel-align-release/accalign -I H -l 32 -t $nthreads -o $output/${aligner}.sam \
$workspace/data/fsva-hg37/hg37.fna $r1 $r2

workspace=/home/yiqing
r1=$workspace/yan/input/BIO250/D3_S1_L001_R1_001.fastq.gz
r2=$workspace/yan/input/BIO250/D3_S1_L001_R2_001.fastq.gz
output=$workspace/yan/output-bio250
aligner=accalign
time $workspace/code/accel-align-release/accalign -I H -l 32 -t $nthreads -o $output/${aligner}.sam \
$workspace/data/fsva-hg37/hg37.fna $r1 $r2