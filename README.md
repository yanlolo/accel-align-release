## Overview ##

Accel-align is a fast alignment tool implemented in C++ programming language.

Download the source code for the various release [here](https://github.com/raja-appuswamy/accel-align-release/releases/)
## Get started ##

### Docker Container ###
You can now pull a preconfigured docker container to get the binaries:
```
docker run -it rajaappuswamy/accel-align
```

### Pre-requirement ###
If you prefer to do a non-docker install, download and install Intel TBB and WFA first. 

#### Intel TBB ####
- source: https://github.com/01org/tbb/releases/tag/2019_U5
- libtbb-dev package

#### WFA ####
- source: https://github.com/smarco/WFA
- Note: please update the `WFA_PATH` in the makefile after install WFA

### Installation ###
* clone it from code repository
* make it: `make -j`

### Build index ###
It's mandatory to build the index before alignment. Options:

```
-l INT the length of k-mers [32]
-m low memory [use if system memory < 16GB]
-t INT number of cpu threads to use[1]
```

Example:

```
path-to-accel-align/accindex -m -l 32 -t 4 path-to-ref/ref.fna
```

It will generate the index aside the reference genome as `path-to-ref/ref.fna.hash`.

### Align ###
When the alignment is triggered, the index will be loaded in memory automatically.

Options:

```
   -t INT number of cpu threads to use[1].
   -l INT length of seed [32].
   -o name of output file to use.
   -x alignment-free.
   -w use WFA for extension. It's using KSW by default.
   -p the maximum distance allowed between the paired-end reads[1000].
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
  