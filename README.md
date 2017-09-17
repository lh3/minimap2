[![Release](https://img.shields.io/badge/Release-v2.1.1-blue.svg?style=flat-square)](https://github.com/lh3/minimap2/releases)
[![BioConda](https://img.shields.io/conda/vn/bioconda/minimap2.svg?style=flat-square)](https://anaconda.org/bioconda/minimap2)
[![PyPI](https://img.shields.io/pypi/v/mappy.svg?style=flat-square)](https://pypi.python.org/pypi/mappy)
[![License](https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square)](LICENSE.txt)
[![Build Status](https://travis-ci.org/lh3/minimap2.svg?branch=master)](https://travis-ci.org/lh3/minimap2)
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
# spliced alignment (no test data)
./minimap2 -ax splice ref.fa rna-seq-reads.fa > spliced.sam
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

If you use minimap2 in your work, please consider to cite:

> Li, H. (2017). Minimap2: fast pairwise alignment for long DNA sequences. [arXiv:1708.01492](https://arxiv.org/abs/1708.01492).

## Installation

For modern x86-64 CPUs, just type `make` in the source code directory. This
will compile a binary `minimap2` which you can copy to your desired location.
If you see compilation errors, try `make sse2only=1` to disable SSE4. Minimap2
will run a little slower. At present, minimap2 does not work with non-x86 CPUs
or ancient CPUs that do not support SSE2. SSE2 is critical to the performance
of minimap2.

## Algorithm Overview

In the following, minimap2 command line options have a dash ahead and are
highlighted in bold.

1. Read **-I** [=*4G*] reference bases, extract (**-k**,**-w**)-minimizers and
   index them in a hash table.

2. Read **-K** [=*200M*] query bases. For each query sequence, do step 3
   through 7:

3. For each (**-k**,**-w**)-minimizer on the query, check against the reference
   index. If a reference minimizer is not among the top **-f** [=*2e-4*] most
   frequent, collect its the occurrences in the reference, which are called
   *seeds*.

4. Sort seeds by position in the reference. Chain them with dynamic
   programming. Each chain represents a potential mapping. For read
   overlapping, report all chains and then go to step 8. For reference mapping,
   do step 5 through 7:

5. Let *P* be the set of primary mappings, which is an empty set initially. For
   each chain from the best to the worst according to their chaining scores: if
   on the query, the chain overlaps with a chain in *P* by **--mask-level**
   [=*0.5*] or higher fraction of the shorter chain, mark the chain as
   *secondary* to the chain in *P*; otherwise, add the chain to *P*.

6. Retain all primary mappings. Also retain up to **-N** [=*5*] top secondary
   mappings if their chaining scores are higher than **-p** [=*0.8*] of their
   corresponding primary mappings.

7. If alignment is requested, filter out an internal seed if it potentially
   leads to both a long insertion and a long deletion. Extend from the
   left-most seed. Perform global alignments between internal seeds.  Split the
   chain if the accumulative score along the global alignment drops by **-z**
   [=*400*], disregarding long gaps. Extend from the right-most seed.  Output
   chains and their alignments.

8. If there are more query sequences in the input, go to step 2 until no more
   queries are left.

9. If there are more reference sequences, reopen the query file from the start
   and go to step 1; otherwise stop.

## Limitations

* Minimap2 may produce suboptimal alignments through long low-complexity
  regions where seed positions may be suboptimal. This should not be a big
  concern because even the optimal alignment may be wrong in such regions.

* Minimap2 does not work well with Illumina short reads as of now.

* Minimap2 requires SSE2 instructions to compile. It is possible to add
  non-SSE2 support, but it would make minimap2 slower by several times.

In general, minimap2 is a young project with most code written since June, 2017.
It may have bugs and room for improvements. Bug reports and suggestions are
warmly welcomed.



[paf]: https://github.com/lh3/miniasm/blob/master/PAF.md
[sam]: https://samtools.github.io/hts-specs/SAMv1.pdf
[minimap]: https://github.com/lh3/minimap
[smartdenovo]: https://github.com/ruanjue/smartdenovo
[longislnd]: https://www.ncbi.nlm.nih.gov/pubmed/27667791
[gaba]: https://github.com/ocxtal/libgaba
[ksw2]: https://github.com/lh3/ksw2
