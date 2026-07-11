CFLAGS=		-g -Wall -O2 -Wc++-compat #-Wextra
CXX=		c++
CXXFLAGS=	-g -Wall -O2 -std=c++17
CPPFLAGS=	-DHAVE_KALLOC
INCLUDES=
OBJS=		kthread.o kalloc.o misc.o bseq.o sketch.o sdust.o options.o index.o \
			lchain.o align.o hit.o seed.o jump.o map.o format.o pe.o esterr.o splitidx.o \
			ksw2_ll_sse.o
PROG=		minimap2
PROG_EXTRA=	sdust minimap2-lite
LIBS=		-lm -lz -lpthread

ifneq ($(aarch64),)
	arm_neon=1
endif

ifeq ($(arm_neon),) # if arm_neon is not defined
ifeq ($(sse2only),) # if sse2only is not defined
	OBJS+=ksw2_extz2_sse41.o ksw2_extd2_sse41.o ksw2_exts2_sse41.o ksw2_extz2_sse2.o ksw2_extd2_sse2.o ksw2_exts2_sse2.o ksw2_dispatch.o
# Wider AVX2 extd2 kernel (mm2-fast port), x86 only. Opt-in: `make avx=1`.
# Runtime dispatch routes the default two-piece extension to AVX2 when available;
# output is byte-identical to SSE (validated). No effect on ARM builds.
# `avx512=1` (implies avx) additionally builds the AVX-512 kernel, which requires
# gcc/icc (clang rejects its non-constant shuffle immediates).
ifneq ($(avx512),)
	avx=1
endif
ifneq ($(avx),)
	OBJS+=ksw2_extd2_avx2.o
	AVX_DISPATCH=-DKSW_AVX_DISPATCH
ifneq ($(avx512),)
	OBJS+=ksw2_extd2_avx512.o
	AVX_DISPATCH+=-DKSW_AVX512_DISPATCH
endif
endif
else                # if sse2only is defined
	OBJS+=ksw2_extz2_sse.o ksw2_extd2_sse.o ksw2_exts2_sse.o
endif
else				# if arm_neon is defined
	OBJS+=ksw2_extz2_neon.o ksw2_extd2_neon.o ksw2_exts2_neon.o
    INCLUDES+=-Isse2neon
ifeq ($(aarch64),)	#if aarch64 is not defined
	CFLAGS+=-D_FILE_OFFSET_BITS=64 -mfpu=neon -fsigned-char
else				#if aarch64 is defined
	CFLAGS+=-D_FILE_OFFSET_BITS=64 -fsigned-char
endif
endif

# Vectorized DP chaining (mm2-fast port). Opt-in: `make vchain=1`.
# On ARM (arm_neon/aarch64) the AVX2 kernel is lowered to NEON via SIMDe;
# on x86 it compiles to native AVX2 (or AVX-512 with avx512=1).
ifneq ($(vchain),)
	OBJS+=mm2fast_chain.o
	CPPFLAGS+=-DMM_VECT_CHAIN
ifneq ($(arm_neon),)   # ARM: use SIMDe to translate AVX2 intrinsics to NEON
	VCHAIN_CXXFLAGS=-DUSE_SIMDE -D__AVX2__ -fsigned-char
	VCHAIN_INCLUDES=-Ilib/simde
else                   # x86: native wide SIMD
ifneq ($(avx512),)
	VCHAIN_CXXFLAGS=-mavx512bw -mavx512f -mavx512dq
else
	VCHAIN_CXXFLAGS=-mavx2
endif
endif
	# C++ runtime needed at link time (final program is linked with $(CC)).
	UNAME_S:=$(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	LIBS+=-lc++
else
	LIBS+=-lstdc++
endif
endif

# Native NEON DP chaining kernel (hand-written 4-wide), ARM only. Opt-in:
# `make aarch64=1 neonchain=1`. Alternative to vchain that avoids SIMDe's 8->4
# AVX2 lowering overhead. Byte-identical; plain C (no C++/SIMDe).
ifneq ($(neonchain),)
	OBJS+=neon_chain.o
	CPPFLAGS+=-DMM_VECT_CHAIN
endif

# Aggressive host-specific tuning for Apple Silicon / other AArch64 CPUs.
# Opt-in: `make aarch64=1 arm_tune=1`. Byte-identical output; ~1.1x end-to-end
# on Apple M-series (the alignment/NEON path benefits from -O3 + native uarch).
ifneq ($(arm_tune),)
	CFLAGS+=-O3 -mcpu=native -flto
	CXXFLAGS+=-O3 -mcpu=native -flto
endif

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -ldl
endif

ifneq ($(tsan),)
	CFLAGS+=-fsanitize=thread
	LIBS+=-fsanitize=thread -ldl
endif

.PHONY:all extra clean depend
.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

extra:all $(PROG_EXTRA)

minimap2:main.o libminimap2.a
		$(CC) $(CFLAGS) main.o -o $@ -L. -lminimap2 $(LIBS)

minimap2-lite:example.o libminimap2.a
		$(CC) $(CFLAGS) $< -o $@ -L. -lminimap2 $(LIBS)

libminimap2.a:$(OBJS)
		$(AR) -csru $@ $(OBJS)

sdust:sdust.c kalloc.o kalloc.h kdq.h kvec.h kseq.h ketopt.h sdust.h
		$(CC) -D_SDUST_MAIN $(CFLAGS) $< kalloc.o -o $@ -lz

# Vectorized chaining wrapper (C++). Built only when vchain=1.
mm2fast_chain.o:mm2fast_chain.cpp contrib/parallel_chaining_v2_22.h mmpriv.h minimap.h
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(VCHAIN_CXXFLAGS) -Icontrib $(VCHAIN_INCLUDES) $(INCLUDES) $< -o $@

# Native NEON chaining kernel (C). Built only when neonchain=1. -ffp-contract=off
# keeps float rounding identical to the reference (no FMA contraction).
neon_chain.o:neon_chain.c mmpriv.h minimap.h
		$(CC) -c $(CFLAGS) -ffp-contract=off $(CPPFLAGS) $(INCLUDES) $< -o $@

# SSE-specific targets on x86/x86_64

ifeq ($(arm_neon),)   # if arm_neon is defined, compile this target with the default setting (i.e. no -msse2)
ksw2_ll_sse.o:ksw2_ll_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse2 $(CPPFLAGS) $(INCLUDES) $< -o $@
endif

ksw2_extz2_sse41.o:ksw2_extz2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

ksw2_extz2_sse2.o:ksw2_extz2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse2 -mno-sse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $(INCLUDES) $< -o $@

ksw2_extd2_sse41.o:ksw2_extd2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

ksw2_extd2_sse2.o:ksw2_extd2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse2 -mno-sse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $(INCLUDES) $< -o $@

ksw2_exts2_sse41.o:ksw2_exts2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

ksw2_exts2_sse2.o:ksw2_exts2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -msse2 -mno-sse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $(INCLUDES) $< -o $@

ksw2_dispatch.o:ksw2_dispatch.c ksw2.h
		$(CC) -c $(CFLAGS) -msse4.1 $(CPPFLAGS) -DKSW_CPU_DISPATCH $(AVX_DISPATCH) $(INCLUDES) $< -o $@

# AVX2 / AVX-512 two-piece extension kernels (x86, avx=1). Same source, ISA-selected
# by __AVX2__ / __AVX512BW__. Byte-identical to SSE (validated).
ksw2_extd2_avx2.o:ksw2_extd2_avx.c ksw2_extd2_avx.h ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -mavx2 $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

ksw2_extd2_avx512.o:ksw2_extd2_avx.c ksw2_extd2_avx.h ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) -mavx512bw -mavx512f -mavx512dq $(CPPFLAGS) -DKSW_CPU_DISPATCH $(INCLUDES) $< -o $@

# NEON-specific targets on ARM

ksw2_extz2_neon.o:ksw2_extz2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_SSE2_ONLY -D__SSE2__ $(INCLUDES) $< -o $@

ksw2_extd2_neon.o:ksw2_extd2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_SSE2_ONLY -D__SSE2__ $(INCLUDES) $< -o $@

ksw2_exts2_neon.o:ksw2_exts2_sse.c ksw2.h kalloc.h
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DKSW_SSE2_ONLY -D__SSE2__ $(INCLUDES) $< -o $@

# other non-file targets

clean:
		rm -fr gmon.out *.o a.out $(PROG) $(PROG_EXTRA) *~ *.a *.dSYM build dist mappy*.so mappy.c python/mappy.c mappy.egg* .eggs

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

# DO NOT DELETE

align.o: minimap.h mmpriv.h bseq.h kseq.h ksw2.h kalloc.h
bseq.o: bseq.h kvec.h kalloc.h kseq.h
esterr.o: mmpriv.h minimap.h bseq.h kseq.h
example.o: minimap.h kseq.h
format.o: kalloc.h mmpriv.h minimap.h bseq.h kseq.h
hit.o: mmpriv.h minimap.h bseq.h kseq.h kalloc.h khash.h
index.o: kthread.h bseq.h minimap.h mmpriv.h kseq.h ksw2.h kalloc.h kvec.h
index.o: khash.h ksort.h
jump.o: mmpriv.h minimap.h bseq.h kseq.h
kalloc.o: kalloc.h
ksw2_extd2_sse.o: ksw2.h kalloc.h
ksw2_exts2_sse.o: ksw2.h kalloc.h
ksw2_extz2_sse.o: ksw2.h kalloc.h
ksw2_ll_sse.o: ksw2.h kalloc.h
kthread.o: kthread.h
lchain.o: mmpriv.h minimap.h bseq.h kseq.h kalloc.h krmq.h
main.o: bseq.h minimap.h mmpriv.h kseq.h ketopt.h
map.o: kthread.h kvec.h kalloc.h sdust.h mmpriv.h minimap.h bseq.h kseq.h
map.o: khash.h ksort.h
misc.o: mmpriv.h minimap.h bseq.h kseq.h ksort.h
options.o: mmpriv.h minimap.h bseq.h kseq.h
pe.o: mmpriv.h minimap.h bseq.h kseq.h kvec.h kalloc.h ksort.h
sdust.o: kalloc.h kdq.h kvec.h sdust.h
seed.o: mmpriv.h minimap.h bseq.h kseq.h kalloc.h ksort.h
sketch.o: kvec.h kalloc.h mmpriv.h minimap.h bseq.h kseq.h
splitidx.o: mmpriv.h minimap.h bseq.h kseq.h
