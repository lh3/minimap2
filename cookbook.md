## Table of Contents

- [Installation](#install)
- [Mapping Genomic Reads](#map-reads)

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

### Map PacBio reads
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
line. *Please always make sure you are using an intended pre-built index.*

### Map Illumina paired-end reads:
```sh
minimap2 -ax sr ecoli_ref.fa ecoli_mason_1.fq ecoli_mason_2.fq > mapped-sr.sam
```

### Evaluating mapping accuracy with simulated reads:
```sh
minimap2 -ax sr ecoli_ref.fa ecoli_mason_1.fq ecoli_mason_2.fq | paftools.js mapeval -
```
The output is:
```
Q       60      19712   0       0.000000000     19712
Q       0       282     219     0.010953286     19994
U       6
```
where a line starting with `Q` gives:
1. Mapping quality (mapQ) threshold
2. Number of mapped reads between this threshold and the previous mapQ threshold.
3. Number of wrong mappings in the same mapQ interval
4. Accumulative mapping error rate
5. Accumulative number of mappings
A `U` line gives the number of unmapped reads (for SAM input only).
