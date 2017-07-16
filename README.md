## Getting Started
```sh
git clone https://gitlab.com/lh3/minimap2
cd minimap2 && make
# long reads against a reference genome
./minimap2 -ax map10k test/MT-human.fa test/MT-orang.fa > test.sam
# create an index first and then map
./minimap2 -x map10k -d MT-human.mmi test/MT-human.fa
./minimap2 -ax map10k MT-human.mmi test/MT-orang.fa > test.sam
# long-read overlap (no test data)
./minimap2 -x ava10k your-reads.fa your-reads.fa > overlaps.paf
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
faster than most other long-read aligners.

Minimap2 is the successor of [minimap][minimap]. It uses a similar
minimizer-based indexing and seeding algorithm, and improves the original
minimap with homopolyer-compressed k-mers (see also [SMARTdenovo][smartdenovo]
and [longISLND][longislnd]), better chaining and the ability to produce CIGAR
with fast extension alignment (see also [libgaba][gaba] and [ksw2][ksw2]).

## Limitations

At the alignment phase, minimap2 performs global alignments between minimizer
hits. If the positions of these minimizer hits are incorrect, the final alignment
may be suboptimal or broken due to the Z-drop heuristics. In addition, in the
event of a long insertion/deletion, minimap2 may split the long event into
a few smaller events. We will address these issues in future.

Minimap2 does not work well with Illumina short reads as of now.

[paf]: https://github.com/lh3/miniasm/blob/master/PAF.md
[sam]: https://samtools.github.io/hts-specs/SAMv1.pdf
[minimap]: https://github.com/lh3/minimap
[smartdenovo]: https://github.com/ruanjue/smartdenovo
[longislnd]: https://www.ncbi.nlm.nih.gov/pubmed/27667791
[gaba]: https://github.com/ocxtal/libgaba
[ksw2]: https://github.com/lh3/ksw2
