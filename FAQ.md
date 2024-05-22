#### 1. Alignment different with option `-a` or `-c`?

Without `-a`, `-c` or `--cs`, minimap2 only finds *approximate* mapping
locations without detailed base alignment. In particular, the start and end
positions of the alignment are imprecise. With one of those options, minimap2
will perform base alignment, which is generally more accurate but is much
slower.

#### 2. How to map Illumina short reads to noisy long reads?

No good solutions. The better approach is to assemble short reads into contigs
and then map noisy reads to contigs.

#### 3. The output SAM doesn't have a header.

By default, minimap2 indexes 4 billion reference bases (4Gb) in a batch and map
all reads against each reference batch. Given a reference longer than 4Gb,
minimap2 is unable to see all the sequences and thus can't produce a correct
SAM header. In this case, minimap2 doesn't output any SAM header. There are two
solutions to this issue. First, you may increase option `-I` to, for example,
`-I8g` to index more reference bases in a batch. This is preferred if your
machine has enough memory. Second, if your machines doesn't have enough memory
to hold the reference index, you can use the `--split-prefix` option in a
command line like:
```sh
minimap2 -ax map-ont --split-prefix=tmp ref.fa reads.fq
```
This second approach uses less memory, but it is slower and requires temporary
disk space.

#### 4. The output SAM is malformatted.

This typically happens when you use nohup to wrap a minimap2 command line.
Nohup is discouraged as it breaks piping. If you have to use nohup, please
specify an output file with option `-o`.

#### 5. How to output one alignment per read?

You can use `--secondary=no` to suppress secondary alignments (aka multiple
mappings), but you can't suppress supplementary alignment (aka split or
chimeric alignment) this way. You can use samtools to filter out these
alignments:
```sh
minimap2 -ax map-out ref.fa reads.fq | samtools view -F0x900
```
However, this is discouraged as supplementary alignment is informative.
