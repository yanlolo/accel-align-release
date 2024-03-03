# echo "sh pe.sh 1>>/home/yiqing/yan/output/pe.log 2>&1"

## align
code_dir=/home/yiqing/code/accel-align-release
in_dir=/home/yiqing/yan/input/simulate/10m-pe-1/
out_dir=/home/yiqing/yan/output
ref=/home/yiqing/data/fsva-hg37/hg37.fna
nthreads=4

### RMI
$code_dir/accalign -l 32 -R -t $nthreads -o $out_dir/rmi.sam $ref \
$in_dir/sv-10m-100-r1.fastq $in_dir/sv-10m-100-r2.fastq
## accuracy check
python $code_dir/utilities/match.py $in_dir/sv-10m-100-pe-align.sam $out_dir/rmi.sam $out_dir/rmi.csv

### BINARY
$code_dir/accalign -l 32 -B -t $nthreads -o $out_dir/binary.sam $ref \
$in_dir/sv-10m-100-r1.fastq $in_dir/sv-10m-100-r2.fastq
## accuracy check
python $code_dir/utilities/match.py $in_dir/sv-10m-100-pe-align.sam $out_dir/binary.sam $out_dir/binary.csv

### HASH
$code_dir/accalign -l 32 -H -t $nthreads -o $out_dir/hash.sam $ref \
$in_dir/sv-10m-100-r1.fastq $in_dir/sv-10m-100-r2.fastq
## accuracy check
python $code_dir/utilities/match.py $in_dir/sv-10m-100-pe-align.sam $out_dir/hash.sam $out_dir/hash.csv
