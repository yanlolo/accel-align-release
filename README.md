## Overview ##

Accel-align is a fast alignment tool implemented in C++ programming language.

## Get started ##

### Docker Container ###

You can now pull a preconfigured docker container to get the binaries:

```
docker run -it rajaappuswamy/accel-align
```

### Pre-requirement ###

If you prefer to do a non-docker install, download and install Intel TBB first.

#### Intel TBB ####

- source: https://github.com/01org/tbb/releases/tag/2019_U5
- libtbb-dev package

### Installation ###

* clone (git clone --recursive https://github.com/raja-appuswamy/accel-align-release)
* Build it: `make`

### Build index ###

It's mandatory to build the index before alignment. Options:

```
-l INT the length of k-mers [32]
```


Example:

```
path-to-accel-align/accindex -l 32 path-to-ref/ref.fna
```

It will generate the index aside the reference genome as `path-to-ref/ref.fna.hash`.

### Align ###

When the alignment is triggered, the index will be loaded in memory automatically.

Options:

```
   -t INT number of cpu threads to use [1].
   -l INT length of seed [32].
   -o name of output file to use.
   -x alignment-free.
   -w use WFA for extension. It's using KSW by default.
   -p the maximum distance allowed between the paired-end reads [1000].
   -d disable embedding, extend all candidates from seeding (this mode is super slow, only for benchmark).
   Note: maximum read length and read name length supported are 512.
```

### Pair-end alignment ###

``` 
path-to-accel-align/accalign options ref.fna read1.fastq read2.fastq
```

Example:

``` 
path-to-accel-align/accalign -l 32 -t 4 -o output-path/out.sam \
path-to-ref/ref.fna input-path/read1.fastq input-path/read2.fastq
``` 

### Single-end alignment ###

``` 
path-to-accel-align/accalign options ref.fna read.fastq
```

Example:

``` 
path-to-accel-align/accalign -l 32 -t 4 -o output-path/out.sam \
path-to-ref/ref.fna input-path/read.fastq
``` 

### Alignment-free  mode ###

Accel-Align does base-to-base align by default. However, Accel-Align supports alignment-free mapping mode where the
position is reported without the CIGAR string.
`-x` option will enable the alignment-free mode.

Example:

``` 
path-to-accel-align/accalign -l 32 -t 4 -x -o output-path/out.sam \
path-to-ref/ref.fna input-path/read.fastq
``` 

### Bisulfite sequencing read alignment (multiple reference) mode ###
#### Index
Index would generate under the same directory of reference as partX, e.g. hg37.fna.hash.part1, hg37.fna.hash.part2.
```
/media/ssd/ngs-data-analysis/code/accel-align-release/accindex \
-l 32 -s /media/ssd/ngs-data-analysis/data/fsva-hg37/hg37.fna 
```

#### Alignment
```
/media/ssd/ngs-data-analysis/code/accel-align-release/accalign \
-l 32 -t 12 -s /media/ssd/ngs-data-analysis/data/fsva-hg37/hg37.fna \
/media/ssd/ngs-data-analysis/yan/input/1m/sv-1m-100-r-GA.fastq > haha.sam
```

#### Test

##### Single-end

##### Pair-end

###### Simulation (simulate PE reads, and mate1 converts C to T, mate2 converts G to A)
```
sed 's/C/T/ig' sv-1m-100-r1.fastq > sv-1m-100-r1-CT.fastq
sed 's/G/A/ig' sv-1m-100-r2.fastq > sv-1m-100-r2-GA.fastq
```

###### Alignment
```
/media/ssd/ngs-data-analysis/code/accel-align-release/accalign \
-l 32 -t 12 -s /media/ssd/ngs-data-analysis/data/fsva-hg37/hg37.fna \
/media/ssd/ngs-data-analysis/yan/input/1m/sv-1m-100-r1-CT.fastq \
/media/ssd/ngs-data-analysis/yan/input/1m/sv-1m-100-r2-GA.fastq \
> /media/ssd/ngs-data-analysis/yan/output/accalign.sam
```



### Citing Accel-align ###
If you use Accel-align in your work, please cite https://doi.org/10.1186/s12859-021-04162-z :

> Yan, Y., Chaturvedi, N. & Appuswamy, R. 
> Accel-Align: a fast sequence mapper and aligner based on the seed–embed–extend method. 
> BMC Bioinformatics 22, 257 (2021).
