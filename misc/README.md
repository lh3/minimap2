## <a name="started"></a>Getting Started

```sh
# install minimap2
git clone https://github.com/lh3/minimap2
cd minimap2 && make
# install the k8 javascript shell
curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
cp k8-0.2.4/k8-`uname -s` k8              # or copy it to a directory on your $PATH
# export PATH="$PATH:`pwd`:`pwd`/misc"    # run this if k8, minimap2 or paftools.js not on your $PATH
minimap2 --cs test/MT-human.fa test/MT-orang.fa | paftools.js view -     # view alignment
minimap2 -c test/MT-human.fa test/MT-orang.fa | paftools.js stat -       # basic alignment statistics
minimap2 -c --cs test/MT-human.fa test/MT-orang.fa \
  | sort -k6,6 -k8,8n | paftools.js call -L15000 -     # calling variants from asm-to-ref alignment
minimap2 -c test/MT-human.fa test/MT-orang.fa \
  | paftools.js liftover -l10000 - <(echo -e "MT_orang\t2000\t5000")     # liftOver
# no test data for the following examples
paftools.js junceval -e anno.gtf splice.sam > out.txt  # compare splice junctions to annotations
paftools.js splice2bed anno.gtf > anno.bed             # convert GTF/GFF3 to BED12
```

## Table of Contents

- [Getting Started](#started)
- [Introduction](#intro)
- [Evaluation](#eval)
  - [Evaluating mapping accuracy with simulated reads](#mapeval)
  - [Evaluating read overlap sensitivity](#oveval)
- [Calling Variants from Assemblies](#asmvar)

## <a name="intro"></a>Introduction

paftools.js is a script that processes alignments in the [PAF format][paf],
such as converting between formats, evaluating mapping accuracy, lifting over
BED files based on alignment, and calling variants from assembly-to-assembly
alignment. This script *requires* the [k8 Javascript shell][k8] to run. On
Linux or Mac, you can download the precompiled k8 binary with:

```sh
curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
cp k8-0.2.4/k8-`uname -s` $HOME/bin/k8  # assuming $HOME/bin in your $PATH
```

It is highly recommended to copy the executable `k8` to a directory on your
`$PATH` such as `/usr/bin/env` can find it. Like python scripts, once you
install `k8`, you can launch paftools.js in one of the two ways:

```sh
path/to/paftools.js             # only if k8 is on your $PATH
k8 path/to/paftools.js
```

In a nutshell, paftools.js has the following commands:

```
Usage: paftools.js <command> [arguments]
Commands:
  view       convert PAF to BLAST-like (for eyeballing) or MAF
  splice2bed convert spliced alignment in PAF/SAM to BED12
  sam2paf    convert SAM to PAF
  delta2paf  convert MUMmer's delta to PAF
  gff2bed    convert GTF/GFF3 to BED12

  stat       collect basic mapping information in PAF/SAM
  liftover   simplistic liftOver
  call       call variants from asm-to-ref alignment with the cs tag
  bedcov     compute the number of bases covered

  mapeval    evaluate mapping accuracy using mason2/PBSIM-simulated FASTQ
  mason2fq   convert mason2-simulated SAM to FASTQ
  pbsim2fq   convert PBSIM-simulated MAF to FASTQ
  junceval   evaluate splice junction consistency with known annotations
  ov-eval    evaluate read overlap sensitivity using read-to-ref mapping
```

paftools.js seamlessly reads both plain text files and gzip'd text files.

## <a name="eval"></a>Evaluation

### <a name="mapeval"></a>Evaluating mapping accuracy with simulated reads

The **pbsim2fq** command of paftools.js converts the MAF output of [pbsim][pbsim]
to FASTQ and encodes the true mapping position in the read name in a format like
`S1_33!chr1!225258409!225267761!-`. Similarly, the **mason2fq** command
converts [mason2][mason2] simulated SAM to FASTQ.

Command **mapeval** evaluates mapped SAM/PAF. Here is example output:

```
Q       60      32478   0       0.000000000     32478
Q       22      16      1       0.000030775     32494
Q       21      43      1       0.000061468     32537
Q       19      73      1       0.000091996     32610
Q       14      66      1       0.000122414     32676
Q       10      27      3       0.000214048     32703
Q       8       14      1       0.000244521     32717
Q       7       13      2       0.000305530     32730
Q       6       46      1       0.000335611     32776
Q       3       10      1       0.000366010     32786
Q       2       20      2       0.000426751     32806
Q       1       248     94      0.003267381     33054
Q       0       31      17      0.003778147     33085
U       3
```

where each Q-line gives the quality threshold, the number of reads mapped with
mapping quality equal to or greater than the threshold, number of wrong
mappings, accumulative mapping error rate and the accumulative number of
mapped reads. The U-line, if present, gives the number of unmapped reads if
they are present in the SAM file.

Suppose the reported mapping coordinate overlap with the true coordinate like
the following:

```
truth:   --------------------
mapper:           ----------------------
         |<- l1 ->|<-- o -->|<-- l2 -->|
```

Let `r=o/(l1+o+l2)`. The reported mapping is considered correct if `r>0.1` by
default.

### <a name="oveval"></a>Evaluating read overlap sensitivity

Command **ov-eval** takes *sorted* read-to-reference alignment and read
overlaps in PAF as input, and evaluates the sensitivity. For example:

```sh
minimap2 -cx map-pb ref.fa reads.fq.gz | sort -k6,6 -k8,8n > reads-to-ref.paf
minimap2 -x ava-pb reads.fq.gz reads.fq.gz > ovlp.paf
k8 ov-eval.js reads-to-ref.paf ovlp.paf
```

## <a name="asmvar"></a>Calling Variants from Haploid Assemblies

The **call** command of paftools.js calls variants from coordinate-sorted
assembly-to-reference alignment. It calls variants from the [cs tag][cs] and
identifies confident/callable regions as those covered by exactly one contig.
Here are example command lines:

```sh
minimap2 -cx asm5 -t8 --cs ref.fa asm.fa > asm.paf  # keeping this file is recommended; --cs required!
sort -k6,6 -k8,8n asm.paf > asm.srt.paf             # sort by reference start coordinate
k8 paftools.js call asm.srt.paf > asm.var.txt
```

Here is sample output:

```
V   chr1    2276040 2276041 1   60  c   g   LJII01000171.1  1217409 1217410 +
V   chr1    2280409 2280410 1   60  a   g   LJII01000171.1  1221778 1221779 +
V   chr1    2280504 2280505 1   60  a   g   LJII01000171.1  1221873 1221874 +
R   chr1    2325140 2436340
V   chr1    2325287 2325287 1   60  -   ct  LJII01000171.1  1272894 1272896 +
V   chr1    2325642 2325644 1   60  tt  -   LJII01000171.1  1273251 1273251 +
V   chr1    2326051 2326052 1   60  c   t   LJII01000171.1  1273658 1273659 +
V   chr1    2326287 2326288 1   60  c   t   LJII01000171.1  1273894 1273895 +
```

where a line starting with `R` gives regions covered by one query contig, and a
V-line encodes a variant in the following format: chr, start, end, query depth,
mapping quality, REF allele, ALT allele, query name, query start, end and the
query orientation. Generally, you should only look at variants where column 5
is one.

By default, when calling variants, "paftools.js call" ignores alignments 50kb
or shorter; when deriving callable regions, it ignores alignments 10kb or
shorter.  It uses two thresholds to avoid edge effects. These defaults are
designed for long-read assemblies. For short reads, both should be reduced.



[paf]: https://github.com/lh3/miniasm/blob/master/PAF.md
[cs]: https://github.com/lh3/minimap2#cs
[k8]: https://github.com/attractivechaos/k8
[maf]: https://genome.ucsc.edu/FAQ/FAQformat#format5
[pbsim]: https://github.com/pfaucon/PBSIM-PacBio-Simulator
[mason2]: https://github.com/seqan/seqan/tree/master/apps/mason2
