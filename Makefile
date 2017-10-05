CFLAGS=		-g -Wall -O2 -Wc++-compat
CPPFLAGS=	-DHAVE_KALLOC
INCLUDES=
OBJS=		kthread.o kalloc.o misc.o bseq.o sketch.o sdust.o index.o chain.o align.o hit.o map.o format.o pe.o ksw2_ll_sse.o
PROG=		minimap2
PROG_EXTRA=	sdust minimap2-lite
LIBS=		-lm -lz -lpthread

ifeq ($(sse2only),)
	OBJS+=ksw2_extz2_sse41.o ksw2_extd2_sse41.o ksw2_exts2_sse41.o ksw2_extz2_sse2.o ksw2_extd2_sse2.o ksw2_exts2_sse2.o ksw2_dispatch.o
else
	OBJS+=ksw2_extz2_sse.o ksw2_extd2_sse.o ksw2_exts2_sse.o
endif

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

extra:all $(PROG_EXTRA)

minimap2:main.o getopt.o libminimap2.a
		$(CC) $(CFLAGS) main.o getopt.o -o $@ -L. -lminimap2 $(LIBS)

minimap2-lite:example.o libminimap2.a
		$(CC) $(CFLAGS) $< -o $@ -L. -lminimap2 $(LIBS)

libminimap2.a:$(OBJS)
		$(AR) -csru $@ $(OBJS)

sdust:sdust.c getopt.o kalloc.o kalloc.h kdq.h kvec.h kseq.h sdust.h
		$(CC) -D_SDUST_MAIN $(CFLAGS) $< getopt.o kalloc.o -o $@ -lz

ksw2_extz2_sse41.o:ksw2_extz2_sse.c ksw2.h kalloc.h
		$(CC) -c -msse4 $(CFLAGS) $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

ksw2_extz2_sse2.o:ksw2_extz2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $(INCLUDES) $< -o $@

ksw2_extd2_sse41.o:ksw2_extd2_sse.c ksw2.h kalloc.h
		$(CC) -c -msse4 $(CFLAGS) $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

ksw2_extd2_sse2.o:ksw2_extd2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $(INCLUDES) $< -o $@

ksw2_exts2_sse41.o:ksw2_exts2_sse.c ksw2.h kalloc.h
		$(CC) -c -msse4 $(CFLAGS) $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

ksw2_exts2_sse2.o:ksw2_exts2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $(INCLUDES) $< -o $@

ksw2_dispatch.o:ksw2_dispatch.c ksw2.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

clean:
		rm -fr gmon.out *.o a.out $(PROG) $(PROG_EXTRA) *~ *.a *.dSYM build dist mappy.so mappy.c python/mappy.c mappy.egg*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

# DO NOT DELETE

align.o: minimap.h mmpriv.h bseq.h ksw2.h kalloc.h
bseq.o: bseq.h kvec.h kalloc.h kseq.h
chain.o: minimap.h mmpriv.h bseq.h kalloc.h
example.o: minimap.h kseq.h
format.o: kalloc.h mmpriv.h minimap.h bseq.h
getopt.o: getopt.h
hit.o: mmpriv.h minimap.h bseq.h kalloc.h khash.h
index.o: kthread.h bseq.h minimap.h mmpriv.h kvec.h kalloc.h khash.h
kalloc.o: kalloc.h
ksw2_extd2_sse.o: ksw2.h kalloc.h
ksw2_exts2_sse.o: ksw2.h kalloc.h
ksw2_extz2_sse.o: ksw2.h kalloc.h
ksw2_ll_sse.o: ksw2.h kalloc.h
main.o: bseq.h minimap.h mmpriv.h getopt.h
map.o: kthread.h kvec.h kalloc.h sdust.h mmpriv.h minimap.h bseq.h khash.h
misc.o: minimap.h ksort.h
pe.o: mmpriv.h minimap.h bseq.h kvec.h kalloc.h ksort.h
sdust.o: kalloc.h kdq.h kvec.h sdust.h
sketch.o: kvec.h kalloc.h minimap.h
