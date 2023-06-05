CC=g++

ifneq ($(DEBUG),)
	CFLAGS=-g -Wall -pthread -O0 -DDBGPRINT -isystem./WFA-paper -L$./WFA-paper/build -std=c++17
else
	CFLAGS=-g -Wall -pthread -O3 -isystem./WFA-paper -mavx2 -L./WFA-paper/build -std=c++17
endif

ACCLDFLAGS=./WFA-paper/build/libwfa.a -lz -ltbb
TARGETS=accindex accalign
CPUSRC=reference.cpp accalign.cpp embedding.cpp ksw2_extz2_sse.c bseq.c index.c kthread.c kalloc.c sketch.c misc.c options.c seed.c
IDXSRC=index.cpp embedding.cpp bseq.c index.c kthread.c kalloc.c sketch.c misc.c options.c indexparameters.cpp io.cpp randstrobes.cpp ./ext/xxhash.c pc.cpp sam.cpp cigar.cpp aln.cpp aligner.cpp fastq.cpp nam.cpp paf.cpp cmdline.cpp readlen.cpp refs.cpp
HEADERS=$(wildcard *.h)
HEADERSHPP=$(wildcard *.hpp)

.PHONY: WFA-paper all
all: WFA-paper ${TARGETS}

WFA-paper:
	$(MAKE) -C WFA-paper clean all

accindex: ${IDXSRC} ${HEADERS} ${HEADERSHPP}
	${CC} -o $@ ${IDXSRC} ${ACCLDFLAGS} ${CFLAGS} 

accalign: ${CPUSRC} ${HEADERS} ${HEADERSHPP}
	${CC} -o $@ ${CPUSRC} ${ACCLDFLAGS} ${CFLAGS}

clean:
	$(MAKE) -C WFA-paper clean
	rm ${TARGETS}
