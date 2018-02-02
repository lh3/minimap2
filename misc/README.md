## <a name="started"></a>Getting Started
```sh
curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
cp k8-0.2.4/k8-`uname -s` k8   # or better copy to a directory on PATH
minimap2 --cs test/MT-*.fa | paf2aln.js - | less   # pretty print base alignment
sam2paf.js aln.sam.gz | less -S                    # convert SAM to PAF
gff2bed.js anno.gtf | less -S                      # convert GTF/GFF3 to BED12
minimap2 -cx splice ref.fa rna-seq.fq | splice2bed.js -   # convert splice aln to BED12
```

## Table of Contents

- [Getting Started](#started)
- [Introduction](#intro)
- [Use Cases](#usage)
  - [Convert PAF to other formats](#paf2aln)
  - [Convert SAM to PAF](#sam2paf)
  - [Convert GTF/GFF3 to BED12 format](#gff2bed)
  - [Convert spliced alignment to BED12](#splice2bed)
  - [Convert spliced alignment to BED12](#eval)

## <a name="intro"></a>Introduction

This directory contains auxiliary scripts for format conversion, mapping
accuracy evaluation and miscellaneous purposes. These scripts *require*
the [k8 Javascript shell][k8] to run. On Linux or Mac, you can download
the precompiled k8 binary with:
```sh
curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
cp k8-0.2.4/k8-`uname -s` k8
```
It is highly recommended to copy the executable `k8` to a directory on your
`PATH` such as `/usr/bin/env` can find it.

## <a name="usage"></a>Use Cases

### <a name="paf2aln"></a>Convert PAF to other formats

Script [paf2aln.js](paf2aln.js) converts PAF with the [cs tag][cs] to
[MAF][maf] or BLAST-like output. It only works with minimap2 output generated
using the `--cs` tag.

### <a name="sam2paf"></a>Convert SAM to PAF

Script [sam2paf.js](sam2paf.js) converts alignments in the SAM format to PAF.

### <a name="gff2bed"></a>Convert GTF/GFF3 to BED12 format

Script [gff2bed.js](gff2bed.js) converts GFF format to 12-column BED format. It
seamlessly works with both GTF and GFF3.

### <a name="splice2bed"></a>Convert spliced alignment to BED12

Script [splice2bed.js](splice2bed.js) converts spliced alignment in SAM or PAF
to 12-column BED format.

### <a name="eval"></a>Evaluating mapping accuracy with simulated reads

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

[cs]: https://github.com/lh3/minimap2#cs
[k8]: https://github.com/attractivechaos/k8
[maf]: https://genome.ucsc.edu/FAQ/FAQformat#format5
[pbsim]: https://github.com/pfaucon/PBSIM-PacBio-Simulator
[mason2]: https://github.com/seqan/seqan/tree/master/apps/mason2
