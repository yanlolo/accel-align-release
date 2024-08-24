CC=g++

ifneq ($(DEBUG),)
	CFLAGS=-g -Wall -pthread -O0 -DDBGPRINT -isystem./WFA-paper -L$./WFA-paper/build -std=c++14
else
	CFLAGS=-g -Wall -pthread -O3 -isystem./WFA-paper -mavx2 -L./WFA-paper/build -std=c++14
endif

ACCLDFLAGS=./WFA-paper/build/libwfa.a -lz -ltbb
TARGETS=accindex accalign
CPUSRC=reference.cpp accalign.cpp embedding.cpp ksw2_extz2_sse.c ./strobe/aligner.cpp ./strobe/cigar.cpp ./strobe/ssw/ssw_cpp.cpp ./strobe/ssw/ssw.c
IDXSRC=index.cpp embedding.cpp
HEADERS=$(wildcard *.h)
HEADERSHPP=$(wildcard *.hpp)

.PHONY: WFA-paper all
all: WFA-paper ${TARGETS}

WFA-paper:
	$(MAKE) -C WFA-paper clean all

accindex: ${IDXSRC} ${HEADERS}
	${CC} -o $@ ${IDXSRC} ${ACCLDFLAGS} ${CFLAGS} 

accalign: ${CPUSRC} ${HEADERS} ${HEADERSHPP}
	${CC} -o $@ ${CPUSRC} ${ACCLDFLAGS} ${CFLAGS}

clean:
	$(MAKE) -C WFA-paper clean
	rm ${TARGETS}
