## <a name="started"></a>Getting Started

```sh
curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
cp k8-0.2.4/k8-`uname -s` $HOME/bin/k8    # assuming $HOME/bin in your $PATH
gff2bed.js anno.gtf | less -S             # convert GTF/GFF3 to BED12 (if k8 installed to $PATH)
k8 gff2bed.js anno.gtf | less -S          # convert GTF/GFF3 to BED12 (if k8 not installed)
sam2paf.js aln.sam.gz | less -S           # convert SAM to PAF
minimap2 --cs test/MT-*.fa | paf2aln.js - | less          # pretty print base alignment
minimap2 -cx splice ref.fa rna-seq.fq | splice2bed.js -   # convert splice aln to BED12
```

## Table of Contents

- [Getting Started](#started)
- [Introduction](#intro)
- [Format Conversion](#conv)
  - [Convert PAF to MAF or BLAST-like format](#paf2aln)
  - [Convert SAM to PAF](#sam2paf)
  - [Convert GTF/GFF3 to BED12](#gff2bed)
  - [Convert spliced alignment to BED12](#splice2bed)
- [Evaluation](#eval)
  - [Evaluating mapping accuracy with simulated reads](#mapeval)
  - [Evaluating read overlap sensitivity](#oveval)
- [Calling Variants from Assemblies](#asmvar)

## <a name="intro"></a>Introduction

This directory contains auxiliary scripts for format conversion, mapping
accuracy evaluation and miscellaneous purposes. These scripts *require*
the [k8 Javascript shell][k8] to run. On Linux or Mac, you can download
the precompiled k8 binary with:

```sh
curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
cp k8-0.2.4/k8-`uname -s` $HOME/bin/k8  # assuming $HOME/bin in your $PATH
```

It is highly recommended to copy the executable `k8` to a directory on your
`$PATH` such as `/usr/bin/env` can find it. Like python or perl scripts, once
you install `k8`, you can launch these k8 scripts either with

```sh
path/to/gff2bed.js anno.gtf.gz
```

or with

```sh
k8 path/to/gff2bed.js anno.gtf
```

All k8 scripts seamlessly work with both plain text files and gzip'd text files.

## <a name="conv"></a>Format Conversion

* <a name="paf2aln"></a>Script [paf2aln.js](paf2aln.js) converts PAF with the
  [cs tag][cs] to [MAF][maf] or BLAST-like output. It only works with minimap2
  output generated using the `--cs` tag.

* <a name="sam2paf"></a>Script [sam2paf.js](sam2paf.js) converts alignments in
  the SAM format to PAF.

* <a name="gff2bed"></a>Script [gff2bed.js](gff2bed.js) converts GFF format to
  12-column BED format. It seamlessly works with both GTF and GFF3.

* <a name="splice2bed"></a>Script [splice2bed.js](splice2bed.js) converts
  spliced alignment in SAM or PAF to 12-column BED format.

## <a name="eval"></a>Evaluation

### <a name="mapeval"></a>Evaluating mapping accuracy with simulated reads

Script [sim-pbsim.js](sim-pbsim.js) converts the MAF output of [pbsim][pbsim]
to FASTQ and encodes the true mapping position in the read name in a format like
`S1_33!chr1!225258409!225267761!-`. Similarly, script
[sim-mason2.js](sim-mason2.js) converts [mason2][mason2] simulated SAM to
FASTQ.

Script [sim-eval.js](sim-eval.js) evaluates mapped SAM/PAF. Here is example output:
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
mapped reads. The U-line gives the number of unmapped reads if they are present
in the SAM file.

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

Script [ov-eval.js](ov-eval.js) takes sorted read-to-reference alignment and
read overlaps in PAF as input, and evaluates the sensitivity. For example:

```sh
minimap2 -cx map-pb ref.fa reads.fq.gz | sort -k6,6 -k8,8n > reads-to-ref.paf
minimap2 -x ava-pb reads.fq.gz reads.fq.gz > ovlp.paf
k8 ov-eval.js reads-to-ref.paf ovlp.paf
```

## <a name="asmvar"></a>Calling Variants from Haploid Assemblies

Command `paftools.js call` calls variants from coordinate-sorted
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
V       chr1    3181702 3181703 1       60      c       t
V       chr1    3181730 3181768 1       60      gtcttacacacggagtcttacacacggtcttacacaca  -
R       chr1    3181796 3260557
V       chr1    3181818 3181822 1       60      tgcg    -
V       chr1    3181831 3181832 1       60      a       g
V       chr1    3181832 3181833 1       60      t       c
V       chr1    3181833 3181834 1       60      t       g
V       chr1    3181874 3181874 1       60      -       ca
V       chr1    3181879 3181880 1       60      g       a
V       chr1    3181886 3181887 1       60      c       g
V       chr1    3181911 3181911 1       60      -       agtcttacacatgcagtcttacacat
V       chr1    3181924 3181925 1       60      t       c
V       chr1    3182079 3182080 1       60      g       a
V       chr1    3182150 3182151 1       60      t       c
V       chr1    3182336 3182337 1       60      t       c
```

where a line starting with `R` gives regions covered by one contig, and a
V-line encodes a variant in the following format: chr, start, end, contig
depth, mapping quality, REF allele and ALT allele.

By default, when calling variants, this script ignores alignments 50kb or
shorter; when deriving callable regions, it ignores alignments 10kb or shorter.
It uses two thresholds to avoid edge effects. These defaults are designed for
long-read assemblies. For short reads, both should be reduced.



[cs]: https://github.com/lh3/minimap2#cs
[k8]: https://github.com/attractivechaos/k8
[maf]: https://genome.ucsc.edu/FAQ/FAQformat#format5
[pbsim]: https://github.com/pfaucon/PBSIM-PacBio-Simulator
[mason2]: https://github.com/seqan/seqan/tree/master/apps/mason2
