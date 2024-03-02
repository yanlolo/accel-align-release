CC=g++

ifneq ($(DEBUG),)
	CFLAGS=-g -Wall -pthread -O0 -DDBGPRINT -isystem./WFA-paper -L$./WFA-paper/build -std=c++17 -I./include
else
	CFLAGS=-g -Wall -pthread -O3 -isystem./WFA-paper -mavx2 -L./WFA-paper/build -std=c++17 -I./include
endif

ACCLDFLAGS=./WFA-paper/build/libwfa.a -lz -ltbb
TARGETS=accindex accalign
CPUSRC=reference.cpp accalign.cpp embedding.cpp ksw2_extz2_sse.c bseq.c index.c kthread.c kalloc.c sketch.c misc.c options.c seed.c strobealign/strobe-index.cpp strobealign/indexparameters.cpp strobealign/io.cpp strobealign/randstrobes.cpp ./strobealign/ext/xxhash.c strobealign/pc.cpp strobealign/sam.cpp strobealign/cigar.cpp strobealign/aln.cpp strobealign/aligner.cpp strobealign/fastq.cpp strobealign/nam.cpp strobealign/paf.cpp strobealign/cmdline.cpp strobealign/readlen.cpp strobealign/refs.cpp ./strobealign/ssw/ssw.c ./strobealign/ssw/ssw_cpp.cpp
IDXSRC=index.cpp embedding.cpp bseq.c index.c kthread.c kalloc.c sketch.c misc.c options.c strobealign/strobe-index.cpp strobealign/indexparameters.cpp strobealign/io.cpp strobealign/randstrobes.cpp ./strobealign/ext/xxhash.c strobealign/pc.cpp strobealign/sam.cpp strobealign/cigar.cpp strobealign/aln.cpp strobealign/aligner.cpp strobealign/fastq.cpp strobealign/nam.cpp strobealign/paf.cpp strobealign/cmdline.cpp strobealign/readlen.cpp strobealign/refs.cpp ./strobealign/ssw/ssw.c ./strobealign/ssw/ssw_cpp.cpp
HEADERS=$(wildcard *.h ./include/*.h)
HEADERSHPP=$(wildcard *.hpp)
RMI_IDXSRC=src/key_gen.cpp

.PHONY: WFA-paper all
all: WFA-paper ${TARGETS}


key_gen: WFA-paper ${RMI_IDXSRC} ${HEADERS}
	${CXX} -o $@ ${RMI_IDXSRC} ${ACCLDFLAGS} ${CFLAGS} -pthread -lstdc++fs

WFA-paper:
	$(MAKE) -C WFA-paper clean all

accindex: ${IDXSRC} ${HEADERS} ${HEADERSHPP}
	${CC} -o $@ ${IDXSRC} ${ACCLDFLAGS} ${CFLAGS} -pthread

accalign: ${CPUSRC} ${HEADERS} ${HEADERSHPP}
	${CC} -o $@ ${CPUSRC} ${ACCLDFLAGS} ${CFLAGS}

clean:
	$(MAKE) -C WFA-paper clean
	rm ${TARGETS}
