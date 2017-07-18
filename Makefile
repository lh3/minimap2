CC=			gcc
CFLAGS=		-g -Wall -O2 -Wc++-compat
CPPFLAGS=	-DHAVE_KALLOC
INCLUDES=	-I.
OBJS=		kthread.o kalloc.o ksw2_extz2_sse.o ksw2_extd2_sse.o misc.o bseq.o \
			sketch.o sdust.o index.o chain.o align.o hit.o map.o format.o
PROG=		minimap2
PROG_EXTRA=	sdust
LIBS=		-lm -lz -lpthread

ifeq ($(sse2only),)
	CFLAGS+=-msse4
endif

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

extra:all $(PROG_EXTRA)

minimap2:main.o libminimap2.a
		$(CC) $(CFLAGS) $< -o $@ -L. -lminimap2 $(LIBS)

minimap2-lite:example.o libminimap2.a
		$(CC) $(CFLAGS) $< -o $@ -L. -lminimap2 $(LIBS)

libminimap2.a:$(OBJS)
		$(AR) -csru $@ $(OBJS)

sdust:sdust.c kalloc.o kalloc.h kdq.h kvec.h kseq.h sdust.h
		$(CC) -D_SDUST_MAIN $(CFLAGS) $< kalloc.o -o $@ -lz

clean:
		rm -fr gmon.out *.o a.out $(PROG) $(PROG_EXTRA) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

# DO NOT DELETE

align.o: minimap.h mmpriv.h bseq.h ksw2.h kalloc.h
bseq.o: bseq.h kseq.h
chain.o: minimap.h mmpriv.h bseq.h kalloc.h
format.o: mmpriv.h minimap.h bseq.h
hit.o: mmpriv.h minimap.h bseq.h kalloc.h
index.o: kthread.h bseq.h minimap.h mmpriv.h kvec.h kalloc.h khash.h
kalloc.o: kalloc.h
ksw2_extd2_sse.o: ksw2.h kalloc.h
ksw2_extz2_sse.o: ksw2.h kalloc.h
main.o: bseq.h minimap.h mmpriv.h
map.o: kthread.h kvec.h kalloc.h sdust.h mmpriv.h minimap.h bseq.h
misc.o: minimap.h ksort.h
sdust.o: kalloc.h kdq.h kvec.h sdust.h
sketch.o: kvec.h kalloc.h minimap.h
