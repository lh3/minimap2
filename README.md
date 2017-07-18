## Getting Started
```sh
git clone https://github.com/lh3/minimap2
cd minimap2 && make
# long reads against a reference genome
./minimap2 -ax map10k test/MT-human.fa test/MT-orang.fa > test.sam
# create an index first and then map
./minimap2 -x map10k -d MT-human.mmi test/MT-human.fa
./minimap2 -ax map10k MT-human.mmi test/MT-orang.fa > test.sam
# long-read overlap (no test data)
./minimap2 -x ava-pb your-reads.fa your-reads.fa > overlaps.paf
# man page
man ./minimap2.1
```

## Introduction

Minimap2 is a fast sequence mapping and alignment program that can find
overlaps between long noisy reads, or map long reads or their assemblies to a
reference genome optionally with detailed alignment (i.e. CIGAR). At present,
it works efficiently with query sequences from a few kilobases to ~100
megabases in length at an error rate ~15%. Minimap2 outputs in the [PAF][paf] or
the [SAM format][sam]. On limited test data sets, minimap2 is over 20 times
faster than most other long-read aligners. It will replace BWA-MEM for long
reads and contig alignment.

Minimap2 is the successor of [minimap][minimap]. It uses a similar
minimizer-based indexing and seeding algorithm, and improves the original
minimap with homopolyer-compressed k-mers (see also [SMARTdenovo][smartdenovo]
and [longISLND][longislnd]), better chaining and the ability to produce CIGAR
with fast extension alignment (see also [libgaba][gaba] and [ksw2][ksw2]) and
piece-wise affine gap cost.

## Installation

For modern x86-64 CPUs, just type `make` in the source code directory. This
will compile a binary `minimap2` which you can copy to your desired location.
If you see compilation errors, try `make sse2only=1` to disable SSE4. Minimap2
will run a little slower. At present, minimap2 does not work with non-x86 CPUs
or ancient CPUs that do not support SSE2. SSE2 is critical to the performance
of minimap2.

## Limitations

* At the alignment phase, minimap2 performs global alignments between minimizer
  hits. If the positions of these minimizer hits are incorrect, the final
  alignment may be suboptimal or unnecessarily fragmented.

* Minimap2 may produce poor alignments that may need post-filtering. We are
  still exploring a reliable and consistent way to report good alignments.

* Minimap2 does not work well with Illumina short reads as of now.

* Minimap2 requires SSE2 instructions to compile. It is possible to add
  non-SSE2 support, but it would make minimap2 slower by several times.

In general, minimap2 is a young project with most code written since June,
2017. It may have bugs and room for improvements. Bug reports and suggestions
are warmly welcomed.



[paf]: https://github.com/lh3/miniasm/blob/master/PAF.md
[sam]: https://samtools.github.io/hts-specs/SAMv1.pdf
[minimap]: https://github.com/lh3/minimap
[smartdenovo]: https://github.com/ruanjue/smartdenovo
[longislnd]: https://www.ncbi.nlm.nih.gov/pubmed/27667791
[gaba]: https://github.com/ocxtal/libgaba
[ksw2]: https://github.com/lh3/ksw2
