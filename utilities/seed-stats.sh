# echo "sh seed-stats.sh 1>>/home/yiqing/yan/output/pipe.log 2>&1"

code_dir=/home/yiqing/code/accel-align-release/
ref=/home/yiqing/data/fsva-hg37/hg37.fna
in_dir=/home/yiqing/yan/input/simulate/1m/
out_dir=/home/yiqing/yan/output/
nthreads=4
rlen=100

#mod=536870911
#xxh=0

for mod in 536870879 536870911 1073741789 1073741823 2147483629 4294967231
do
  for xxh in 0 32 64
  do
    ## index stats
    mod=4294967231
    xxh=0
    $code_dir/accindex_stats -l 32 -h ${mod} -x ${xxh} $ref
    python $code_dir/utilities/seed-stats.py ${mod}-${xxh} $out_dir

    $code_dir/accindex -l 32 -h ${mod} -x ${xxh} $ref

    time $code_dir/accalign -l 32 -H -t $nthreads -o $out_dir/hash.sam $ref \
    $in_dir/sv-10m-${rlen}-r1.fastq $in_dir/sv-10m-${rlen}-r2.fastq

    ## accuracy check
    python $code_dir/utilities/match.py $in_dir/sv-1m-${rlen}-pe-align.sam $out_dir/hash.sam $out_dir/hash.csv
  done
done