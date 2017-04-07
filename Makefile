CC=			gcc
CFLAGS=		-g -Wall -O2 -Wc++-compat -Wno-unused-function
CPPFLAGS=
INCLUDES=	-I.
OBJS=		kalloc.o kthread.o bseq.o sketch.o sdust.o index.o
PROG=		minimap2
PROG_EXTRA=	sdust
LIBS=		-lm -lz -lpthread

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

extra:all $(PROG_EXTRA)

minimap2:libminimap2.a
		$(CC) $(CFLAGS) $< -o $@ -L. -lminimap2 $(LIBS)

libminimap2.a:$(OBJS)
		$(AR) -csru $@ $(OBJS)

sdust:sdust.c kalloc.o kdq.h kvec.h kseq.h sdust.h
		$(CC) -D_SDUST_MAIN $(CFLAGS) $< kalloc.o -o $@ -lz

clean:
		rm -fr gmon.out *.o a.out $(PROG) $(PROG_EXTRA) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

bseq.o: bseq.h kseq.h
index.o: kthread.h bseq.h minimap.h kvec.h kalloc.h khash.h
kalloc.o: kalloc.h
misc.o: minimap.h ksort.h
sdust.o: kalloc.h kdq.h kvec.h sdust.h
sketch.o: kvec.h kalloc.h minimap.h
