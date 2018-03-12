## Table of Contents

- [Installation](#install)
- [Mapping Genomic Reads](#map-reads)

## <a name="install"></a>Installation

```sh
curl -L https://github.com/lh3/minimap2/releases/download/v2.9/minimap2-2.9_x64-linux.tar.bz2 | tar jxf -
export PATH="$PATH:"`pwd`                              # set up PATH
cp minimap2-2.9_x64-linux/{minimap2,k8,paftools.js} .  # copy executables
cp minimap2-2.9_x64-linux/test/MT-*.fa .               # copy small examples
curl -L https://github.com/lh3/minimap2/releases/download/v2.0/ecoli.tgz | tar zxf -
```

## <a name="map-reads"></a>Mapping Genomic Reads

* Map example E. coli reads (takes about 12 wall-clock seconds):
  ```sh
  minimap2 -ax map-pb -t4 ecoli_ref.fa ecoli_p6_25x_canu.fa > mapped.sam
  ```
