# echo "sh pe.sh 1>>/home/yiqing/yan/output/pe.log 2>&1"

## align
code_dir=/home/yiqing/code/
in_dir=/home/yiqing/yan/input/simulate/1m/
out_dir=/home/yiqing/yan/output
ref=/home/yiqing/data/fsva-hg37/hg37.fna
nthreads=12
rlen=100
r=sv-1m-${rlen}-r.fastq
truth=sv-1m-${rlen}-se-align.sam

########### Minimap2 ###########
aligner=minimap2
time $code_dir/minimap2/minimap2 -t $nthreads -ax sr /home/yiqing/data/fsva-hg37/hg37.mmi \
$in_dir/$r > $out_dir/${aligner}.sam
python $code_dir/accel-align-release/utilities/match.py $in_dir/$truth $out_dir/${aligner}.sam $out_dir/${aligner}.csv

########### strobeAlign ###########
aligner=strobealign
time $code_dir/strobealign/build/strobealign --use-index $ref -t $nthreads $in_dir/$r  > $out_dir/${aligner}.sam
python $code_dir/accel-align-release/utilities/match.py $in_dir/$truth $out_dir/${aligner}.sam $out_dir/${aligner}.csv

########### bwa ###########
aligner=bwa
time $code_dir/bwa/bwa mem -M -t $nthreads $ref $in_dir/$r  > $out_dir/${aligner}.sam
python $code_dir/accel-align-release/utilities/match.py $in_dir/$truth $out_dir/${aligner}.sam $out_dir/${aligner}.csv

########### AA master branch ###########
aligner=accalign
time $code_dir/accel-align-release/accalign -l 32 -t $nthreads -o $out_dir/${aligner}.sam $ref $in_dir/$r
python $code_dir/accel-align-release/utilities/match.py $in_dir/$truth $out_dir/${aligner}.sam $out_dir/${aligner}.csv

########### AA benmark-all branch ###########
aligner=accalign
time $code_dir/accel-align-release/accalign -I H -l 32 -t $nthreads -o $out_dir/${aligner}.sam $ref $in_dir/$r
python $code_dir/accel-align-release/utilities/match.py $in_dir/$truth $out_dir/${aligner}.sam $out_dir/${aligner}.csv
