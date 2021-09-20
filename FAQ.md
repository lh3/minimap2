#### 1. Alignment different with option `-a` or `-c`?

Without `-a`, `-c` or `--cs`, minimap2 only finds *approximate* mapping
locations without detailed base alignment. In particular, the start and end
positions of the alignment are impricise. With one of those options, minimap2
will perform base alignment, which is generally more accurate but is much
slower.

#### 2. How to map Illumina short reads to noisy long reads?

No good solutions. The better approach is to assemble short reads into contigs
and then map noisy reads to contigs.

#### 3. The output SAM doesn't have a header.

By default, minimap2 indexes 4 billion reference bases (4Gb) in a batch and map
all reads against each reference batch. Given a reference longer than 4Gb,
minimap2 is unable to see all the sequences and thus can't produce a correct
SAM header. In this case, minimap2 doesn't output any SAM header. You may
increase option `-I` to, for example, `-I8g` to index more reference bases in a
batch. This is preferred if your machine has enough memory. Otherwise, see the
next question.

#### 4. The reference index is too large to fit in memory.

When the reference is longer than 4Gb (or as set by `-I` described above),
minimap2 processes it in batches and reports mappings for each batch
consecutively. The resulting output stream may therefore seem to report
multiple primary mappings for many reads.

The `--split-prefix` option adds a postprocessing stage to merge results across
the batches and produce mappings more similar (but not exactly identical) to
those for one large reference index.

```sh
minimap2 -ax map-ont --split-prefix=tmp ref.fa reads.fq
```

This writes out temporary files prefixed by `tmp`, using time and disk space.
See [Gamaarachchi et al. (2019)](https://www.nature.com/articles/s41598-019-40739-8)
for further information about this feature.

**Distributed mapping.** The postprocessing merge algorithm can also be used on
intermediate mappings for separate index shards, allowing parallelization on a
cluster. To do this, first generate the reference index shards using
`minimap2 -ax sr -I 9999G -d N.idx ref.N.fa` for each *N* of the desired
shards. Then, map the reads against each shard using
`minimap2 -ax sr --split-map N.mappings N.idx reads_1.fq reads_2.fq` for each
*N* (in parallel as desired). Lastly, merge the intermediate mappings into a
final PAF or SAM with
`minimap2 -ax sr --split-merge reads_1.fq reads_2.fq . *.mappings`. Make sure
to supply the same configuration options to each invocation.

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
