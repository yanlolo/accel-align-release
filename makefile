CC=g++

ifneq ($(DEBUG),)
	CFLAGS=-g -Wall -pthread -O0 -DDBGPRINT -isystem./WFA-paper -L$./WFA-paper/build -std=c++17 -I./include
else
	CFLAGS=-g -Wall -pthread -O3 -isystem./WFA-paper -mavx2 -L./WFA-paper/build -std=c++17 -I./include
endif

ACCLDFLAGS=./WFA-paper/build/libwfa.a -lz -ltbb
TARGETS=accindex accalign
CPUSRC=src/rmi.cpp src/reference.cpp src/accalign.cpp src/embedding.cpp src/ksw2_extz2_sse.c src/bseq.c src/index.c src/kthread.c src/kalloc.c src/sketch.c src/misc.c src/options.c src/seed.c strobealign/strobe-index.cpp strobealign/indexparameters.cpp strobealign/io.cpp strobealign/randstrobes.cpp ./strobealign/ext/xxhash.c strobealign/pc.cpp strobealign/sam.cpp strobealign/cigar.cpp strobealign/aln.cpp strobealign/aligner.cpp strobealign/fastq.cpp strobealign/nam.cpp strobealign/paf.cpp strobealign/cmdline.cpp strobealign/readlen.cpp strobealign/refs.cpp ./strobealign/ssw/ssw.c ./strobealign/ssw/ssw_cpp.cpp
IDXSRC=src/index.cpp src/embedding.cpp src/bseq.c src/index.c src/kthread.c src/kalloc.c src/sketch.c src/misc.c src/options.c strobealign/strobe-index.cpp strobealign/indexparameters.cpp strobealign/io.cpp strobealign/randstrobes.cpp ./strobealign/ext/xxhash.c strobealign/pc.cpp strobealign/sam.cpp strobealign/cigar.cpp strobealign/aln.cpp strobealign/aligner.cpp strobealign/fastq.cpp strobealign/nam.cpp strobealign/paf.cpp strobealign/cmdline.cpp strobealign/readlen.cpp strobealign/refs.cpp ./strobealign/ssw/ssw.c ./strobealign/ssw/ssw_cpp.cpp
STATSIDXSRC=src/index-stats.cpp src/embedding.cpp src/bseq.c src/index.c src/kthread.c src/kalloc.c src/sketch.c src/misc.c src/options.c strobealign/strobe-index.cpp strobealign/indexparameters.cpp strobealign/io.cpp strobealign/randstrobes.cpp ./strobealign/ext/xxhash.c strobealign/pc.cpp strobealign/sam.cpp strobealign/cigar.cpp strobealign/aln.cpp strobealign/aligner.cpp strobealign/fastq.cpp strobealign/nam.cpp strobealign/paf.cpp strobealign/cmdline.cpp strobealign/readlen.cpp strobealign/refs.cpp ./strobealign/ssw/ssw.c ./strobealign/ssw/ssw_cpp.cpp
HEADERS=$(wildcard *.h ./include/*.h)
HEADERSHPP=$(wildcard *.hpp)
RMI_IDXSRC=src/key_gen.cpp
STATSSRC=src/stats.cpp

.PHONY: WFA-paper all
all: WFA-paper ${TARGETS}


key_gen: WFA-paper ${RMI_IDXSRC} ${HEADERS}
	${CXX} -o $@ ${RMI_IDXSRC} ${ACCLDFLAGS} ${CFLAGS} -pthread -lstdc++fs

WFA-paper:
	$(MAKE) -C WFA-paper clean all

accindex: ${IDXSRC} ${HEADERS} ${HEADERSHPP}
	${CC} -o $@ ${IDXSRC} ${ACCLDFLAGS} ${CFLAGS} -pthread -lstdc++fs

accindex_stats: WFA-paper ${STATSIDXSRC} ${HEADERS} ${HEADERSHPP}
	${CC} -o $@ ${STATSIDXSRC} ${ACCLDFLAGS} ${CFLAGS} -pthread

accalign: ${CPUSRC} ${HEADERS} ${HEADERSHPP}
	${CC} -o $@ ${CPUSRC} ${ACCLDFLAGS} ${CFLAGS}

stats: WFA-paper ${STATSSRC} ${HEADERS}
	${CXX} -o $@ ${STATSSRC} ${ACCLDFLAGS} ${CFLAGS} -fopenmp

clean:
	$(MAKE) -C WFA-paper clean
	rm ${TARGETS}
