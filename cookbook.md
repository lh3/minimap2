## Table of Contents

- [Installation](#install)
- [Mapping Genomic Reads](#map-reads)
  * [Mapping long reads](#map-pb)
  * [Mapping Illumina paired-end reads](#map-sr)
  * [Evaluating mapping accuracy with simulated reads (for developers)](#mapeval)
- [Full-Genome Alignment](#genome-aln)
  * [Intra-species assembly alignment](#asm-to-ref)
  * [Calling variants from assembly-to-reference alignment](#asm-var)
- [Read Overlap](#read-overlap)
  * [Long-read overlap](#long-read-overlap)
  * [Evaluating overlap sensitivity (for developers)](#ov-eval)

## <a name="install"></a>Installation

```sh
# install minimap2 executables
curl -L https://github.com/lh3/minimap2/releases/download/v2.9/minimap2-2.9_x64-linux.tar.bz2 | tar jxf -
cp minimap2-2.9_x64-linux/{minimap2,k8,paftools.js} .  # copy executables
export PATH="$PATH:"`pwd`                              # put the current directory on PATH
# download example datasets
curl -L https://github.com/lh3/minimap2/releases/download/v2.0/cookbook-data.tgz | tar zxf -
```

## <a name="map-reads"></a>Mapping Genomic Reads

### <a name="map-pb"></a>Mapping long reads
```sh
minimap2 -ax map-pb -t4 ecoli_ref.fa ecoli_p6_25x_canu.fa > mapped.sam
```
Alternatively, you can create a minimap2 index first and then map:
```sh
minimap2 -x map-pb -d ecoli-pb.mmi ecoli_ref.fa                      # create an index
minimap2 -ax map-pb ecoli-pb.mmi ecoli_p6_25x_canu.fa > mapped.sam
```
This will save you a couple of minutes when you map against the human genome.
**HOWEVER**, key algorithm parameters such as the k-mer length and window
size can't be changed after indexing. Minimap2 will give you a warning if
parameters used in a pre-built index doesn't match parameters on the command
line. **Please always make sure you are using an intended pre-built index.**

### <a name="map-sr"></a>Mapping Illumina paired-end reads:
```sh
minimap2 -ax sr -t4 ecoli_ref.fa ecoli_mason_1.fq ecoli_mason_2.fq > mapped-sr.sam
```

### <a name="mapeval"></a>Evaluating mapping accuracy with simulated reads (for developers)
```sh
minimap2 -ax sr ecoli_ref.fa ecoli_mason_1.fq ecoli_mason_2.fq | paftools.js mapeval -
```
The output is:
```
Q       60      19712   0       0.000000000     19712
Q       0       282     219     0.010953286     19994
U       6
```
where a `U`-line gives the number of unmapped reads (for SAM input only); a
`Q`-line gives:

1. Mapping quality (mapQ) threshold
2. Number of mapped reads between this threshold and the previous mapQ threshold.
3. Number of wrong mappings in the same mapQ interval
4. Accumulative mapping error rate
5. Accumulative number of mappings

For `paftools.js mapeval` to work, you need to encode the true read positions
in read names in the right format. For [PBSIM][pbsim] and [mason2][mason2], we
provide scripts to generate the right format. Simulated reads in this cookbook
were created with the following command lines:
```sh
# in PBSIM source code directory:
src/pbsim ../ecoli_ref.fa --depth 1 --sample-fastq sample/sample.fastq
paftools.js pbsim2fq ../ecoli_ref.fa.fai sd_0001.maf > ../ecoli_pbsim.fa

# mason2 simulation
mason_simulator --illumina-prob-mismatch-scale 2.5 -ir ecoli_ref.fa -n 10000 -o tmp-l.fq -or tmp-r.fq -oa tmp.sam
paftools.js mason2fq tmp.sam | seqtk seq -1 > ecoli_mason_1.fq
paftools.js mason2fq tmp.sam | seqtk seq -1 > ecoli_mason_2.fq
```



## <a name="genome-aln"></a>Full-Genome Alignment

### <a name="asm-to-ref"></a>Intra-species assembly alignment
```sh
minimap2 -cx asm5 --cs ecoli_ref.fa ecoli_canu.fa > ecoli_canu.paf
```
Here `ecoli_canu.fa` is the Canu assembly of `ecoli_p6_25x_canu.fa`. This
command line outputs alignments in the [PAF format][paf]. Use `-a` instead of
`-c` to get output in the SAM format.

### <a name="asm-var"></a>Calling variants from assembly-to-reference alignment
```sh
# don't forget the "--cs" option; otherwise it doesn't work
minimap2 -cx asm5 --cs ecoli_ref.fa ecoli_canu.fa \
  | sort -k6,6 -k8,8n \
  | paftools.js call - > out.txt
```
The first few lines in `out.txt` are
```
R       NC_000913       0       257899
V       NC_000913       668     669     1       60      g       -       tig00000001     3439623 3439623 +
V       NC_000913       896     897     1       60      g       -       tig00000001     3439850 3439850 +
V       NC_000913       1304    1305    1       60      g       -       tig00000001     3440257 3440257 +
```
where an `R`-line gives the regions coverged by one contig; a `V`-line consists of:

1. Reference sequence name
2. Reference start
3. Reference end (a 0-based, half-close-half-open interval)
4. Contig coverage (please ignore a `V`-line if this number is not one!)
5. Mean mapping quality of contigs covering the site
6. Reference base ("-" for an insertion)
7. Contig base ("-" for a deletion)
8. Contig name
9. Contig start
10. Contig end
11. Contig strand

`paftools.js call` does not support VCF output. This feature is planned but has
not been implemented yet.




## <a name="read-overlap"></a>Read Overlap

### <a name="long-read-overlap"></a>Long read overlap
```sh
# For pacbio reads:
minimap2 -x ava-pb ecoli_p6_25x_canu.fa ecoli_p6_25x_canu.fa > overlap.paf
# For Nanopore reads (ava-ont also works with PacBio but not as good):
minimap2 -x ava-ont -r 10000 ecoli_p6_25x_canu.fa ecoli_p6_25x_canu.fa > overlap.paf
# If you have miniasm installed:
miniasm -f ecoli_p6_25x_canu.fa overlap.paf > asm.gfa
```
Here we explicitly applied `-r 10000`. We are considering to set this as the
default for the `ava-ont` mode as this seems to improve the contiguity for
nanopore read assembly (Loman, personal communication).

*Minimap2 doesn't work well with short-read overlap.*

### <a name="ov-eval"></a>Evaluating overlap sensitivity (for developers)

```sh
# read to reference mapping
minimap2 -cx map-pb ecoli_ref.fa ecoli_p6_25x_canu.fa > to-ref.paf
# evaluate overlap sensitivity
sort -k6,6 -k8,8n to-ref.paf | paftools.js ov-eval - overlap.paf
```
You can see that for PacBio reads, minimap2 achieves higher overlap sensitivity
with `-x ava-pb` (99% vs 93% with `-x ava-ont`).



## <a name="map-rna"></a>Mapping Long RNA-seq Reads

(Unfinished section. Data to come later...)

### <a name="map-direct-rna"></a>Mapping Nanopore direct-RNA reads

```sh
minimap2 -ax splice -k14 -uf ref.fa reads.fa > aln.sam
```



[pbsim]: https://github.com/pfaucon/PBSIM-PacBio-Simulator
[mason2]: https://github.com/seqan/seqan/tree/master/apps/mason2
[paf]: https://github.com/lh3/miniasm/blob/master/PAF.md
