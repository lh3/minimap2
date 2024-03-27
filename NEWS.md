Release 2.28-r1209 (27 March 2024)
----------------------------------

Notable changes to minimap2:

 * Bugfix: `--MD` was not working properly due to the addition of `--ds` in the
   last release (#1181 and #1182).

 * New feature: added an experimental preset `lq:hqae` for aligning accurate
   long reads back to their assembly. It has been observed that `map-hifi` and
   `lr:hq` may produce many wrong alignments around centromeres when accurate
   long reads (PacBio HiFi or Nanopore duplex/Q20+) are mapped to a diploid
   assembly constructed from them. This new preset produces much more accurate
   alignment. It is still experimental and may be subjective to changes in
   future.

 * Change: reduced the default `--cap-kalloc` to 500m to lower the peak
   memory consumption (#855).

Notable changes to mappy:

 * Bugfix: mappy option struct was out of sync with minimap2 (#1177).

Minimap2 should output identical alignments to v2.27.

(2.28: 27 March 2024, r1209)



Release 2.27-r1193 (12 March 2024)
----------------------------------

Notable changes to minimap2:

 * New feature: added the `lr:hq` preset for accurate long reads at ~1% error
   rate. This was suggested by Oxford Nanopore developers (#1127). It is not
   clear if this preset also works well for PacBio HiFi reads.

 * New feature: added the `map-iclr` preset for Illumina Complete Long Reads
   (#1069), provided by Illumina developers.

 * New feature: added option `-b` to specify mismatch penalty for base
   transitions (i.e. A-to-G or C-to-T changes).

 * New feature: added option `--ds` to generate a new `ds:Z` tag that
   indicates uncertainty in INDEL positions. It is an extension to `cs`. The
   `mgutils-es6.js` script in minigraph parses `ds`.

 * Bugfix: avoided a NULL pointer dereference (#1154). This would not have an
   effect on most systems but would still be good to fix.

 * Bugfix: reverted the value of `ms:i` to pre-2.22 versions (#1146). This was
   an oversight. See fcd4df2 for details.

Notable changes to paftools.js and mappy:

 * New feature: expose `bw_long` to mappy's Aligner class (#1124).

 * Bugfix: fixed several compatibility issues with k8 v1.0 (#1161 and #1166).
   Subcommands "call", "pbsim2fq" and "mason2fq" were not working with v1.0.

Minimap2 should output identical alignments to v2.26, except the ms tag.

(2.27: 12 March 2024, r1193)



Release 2.26-r1175 (29 April 2023)
----------------------------------

Fixed the broken Python package. This is the only change.

(2.26: 25 April 2023, r1173)



Release 2.25-r1173 (25 April 2023)
----------------------------------

Notable changes:

 * Improvement: use the miniprot splice model for RNA-seq alignment by default.
   This model considers non-GT-AG splice sites and leads to slightly higher
   (<0.1%) accuracy and sensitivity on real human data.

 * Change: increased the default `-I` to `8G` such that minimap2 would create a
   uni-part index for a pair of mammalian genomes. This change may increase the
   memory for all-vs-all read overlap alignment given large datasets.

 * New feature: output the sequences in secondary alignments with option
   `--secondary-seq` (#687).

 * Bugfix: --rmq was not parsed correctly (#1010)

 * Bugfix: possibly incorrect coordinate when applying end bonus to the target
   sequence (#1025). This is a ksw2 bug. It does not affect minimap2 as
   minimap2 is not using the affected feature.

 * Improvement: incorporated several changes for better compatibility with
   Windows (#1051) and for minimap2 integration at Oxford Nanopore Technologies
   (#1048 and #1033).

 * Improvement: output the HD-line in SAM output (#1019).

 * Improvement: check minimap2 index file in mappy to prevent segmentation
   fault for certain indices (#1008).

For genomic sequences, minimap2 should give identical output to v2.24.
Long-read RNA-seq alignment may occasionally differ from previous versions.

(2.25: 25 April 2023, r1173)



Release 2.24-r1122 (26 December 2021)
-------------------------------------

This release improves alignment around long poorly aligned regions. Older
minimap2 may chain through such regions in rare cases which may result in
missing alignments later. The issue has become worse since the the change of
the chaining algorithm in v2.19. v2.23 implements an incomplete remedy. This
release provides a better solution with a X-drop-like heuristic and by enabling
two-bandwidth chaining in the assembly mode.

(2.24: 26 December 2021, r1122)



Release 2.23-r1111 (18 November 2021)
-------------------------------------

Notable changes:

 * Bugfix: fixed missing alignments around long inversions (#806 and #816).
   This bug affected v2.19 through v2.22.

 * Improvement: avoid extremely long mapping time for pathologic reads with
   highly repeated k-mers not in the reference (#771). Use --q-occ-frac=0
   to disable the new heuristic.

 * Change: use --cap-kalloc=1g by default.

(2.23: 18 November 2021, r1111)



Release 2.22-r1101 (7 August 2021)
----------------------------------

When choosing the best alignment, this release uses logarithm gap penalty and
query-specific mismatch penalty. It improves the sensitivity to long INDELs in
repetitive regions.

Other notable changes:

 * Bugfix: fixed an indirect memory leak that may waste a large amount of
   memory given highly repetitive reference such as a 16S RNA database (#749).
   All versions of minimap2 have this issue.

 * New feature: added --cap-kalloc to reduce the peak memory. This option is
   not enabled by default but may become the default in future releases.

Known issue:

 * Minimap2 may take a long time to map a read (#771). So far it is not clear
   if this happens to v2.18 and earlier versions.

(2.22: 7 August 2021, r1101)



Release 2.21-r1071 (6 July 2021)
--------------------------------

This release fixed a regression in short-read mapping introduced in v2.19
(#776). It also fixed invalid comparisons of uninitialized variables, though
these are harmless (#752). Long-read alignment should be identical to v2.20.

(2.21: 6 July 2021, r1071)



Release 2.20-r1061 (27 May 2021)
--------------------------------

This release fixed a bug in the Python module and improves the command-line
compatibiliity with v2.18. In v2.19, if `-r` is specified with an `asm*` preset,
users would get alignments more fragmented than v2.18. This could be an issue
for existing pipelines specifying `-r`. This release resolves this issue.

(2.20: 27 May 2021, r1061)



Release 2.19-r1057 (26 May 2021)
--------------------------------

This release includes a few important improvements backported from unimap:

 * Improvement: more contiguous alignment through long INDELs. This is enabled
   by the minigraph chaining algorithm. All `asm*` presets now use the new
   algorithm. They can find INDELs up to 100kb and may be faster for
   chromosome-long contigs. The default mode and `map*` presets use this
   algorithm to replace the long-join heuristic.

 * Improvement: better alignment in highly repetitive regions by rescuing
   high-occurrence seeds. If the distance between two adjacent seeds is too
   large, attempt to choose a fraction of high-occurrence seeds in-between.
   Minimap2 now produces fewer clippings and alignment break points in long
   satellite regions.

 * Improvement: allow to specify an interval of k-mer occurrences with `-U`.
   For repeat-rich genomes, the automatic k-mer occurrence threshold determined
   by `-f` may be too large and makes alignment impractically slow. The new
   option protects against such cases. Enabled for `asm*` and `map-hifi`.

 * New feature: added the `map-hifi` preset for maping PacBio High-Fidelity
   (HiFi) reads.

 * Change to the default: apply `--cap-sw-mem=100m` for genomic alignment.

 * Bugfix: minimap2 could not generate an index file with `-xsr` (#734).

This release represents the most signficant algorithmic change since v2.1 in
2017. With features backported from unimap, minimap2 now has similar power to
unimap for contig alignment. Unimap will remain an experimental project and is
no longer recommended over minimap2. Sorry for reverting the recommendation in
short time.

(2.19: 26 May 2021, r1057)



Release 2.18-r1015 (9 April 2021)
---------------------------------

This release fixes multiple rare bugs in minimap2 and adds additional
functionality to paftools.js.

Changes to minimap2:

 * Bugfix: a rare segfault caused by an off-by-one error (#489)

 * Bugfix: minimap2 segfaulted due to an uninitilized variable (#622 and #625).

 * Bugfix: minimap2 parsed spaces as field separators in BED (#721). This led
   to issues when the BED name column contains spaces.

 * Bugfix: minimap2 `--split-prefix` did not work with long reference names
   (#394).

 * Bugfix: option `--junc-bonus` didn't work (#513)

 * Bugfix: minimap2 didn't return 1 on I/O errors (#532)

 * Bugfix: the `de:f` tag (sequence divergence) could be negative if there were
   ambiguous bases

 * Bugfix: fixed two undefined behaviors caused by calling memcpy() on
   zero-length blocks (#443)

 * Bugfix: there were duplicated SAM @SQ lines if option `--split-prefix` is in
   use (#400 and #527)

 * Bugfix: option -K had to be smaller than 2 billion (#491). This was caused
   by a 32-bit integer overflow.

 * Improvement: optionally compile against SIMDe (#597). Minimap2 should work
   with IBM POWER CPUs, though this has not been tested. To compile with SIMDe,
   please use `make -f Makefile.simde`.

 * Improvement: more informative error message for I/O errors (#454) and for
   FASTQ parsing errors (#510)

 * Improvement: abort given malformatted RG line (#541)

 * Improvement: better formula to estimate the `dv:f` tag (approximate sequence
   divergence). See DOI:10.1101/2021.01.15.426881.

 * New feature: added the `--mask-len` option to fine control the removal of
   redundant hits (#659). The default behavior is unchanged.

Changes to mappy:

 * Bugfix: mappy caused segmentation fault if the reference index is not
   present (#413).

 * Bugfix: fixed a memory leak via 238b6bb3

 * Change: always require Cython to compile the mappy module (#723). Older
   mappy packages at PyPI bundled the C source code generated by Cython such
   that end users did not need to install Cython to compile mappy. However, as
   Python 3.9 is breaking backward compatibility, older mappy does not work
   with Python 3.9 anymore. We have to add this Cython dependency as a
   workaround.

Changes to paftools.js:

 * Bugfix: the "part10-" line from asmgene was wrong (#581)

 * Improvement: compatibility with GTF files from GenBank (#422)

 * New feature: asmgene also checks missing multi-copy genes

 * New feature: added the misjoin command to evaluate large-scale misjoins and
   megabase-long inversions.

Although given the many bug fixes and minor improvements, the core algorithm
stays the same. This version of minimap2 produces nearly identical alignments
to v2.17 except very rare corner cases.

Now unimap is recommended over minimap2 for aligning long contigs against a
reference genome. It often takes less wall-clock time and is much more
sensitive to long insertions and deletions.

(2.18: 9 April 2021, r1015)



Release 2.17-r941 (4 May 2019)
------------------------------

Changes since the last release:

 * Fixed flawed CIGARs like `5I6D7I` (#392).

 * Bugfix: TLEN should be 0 when either end is unmapped (#373 and #365).

 * Bugfix: mappy is unable to write index (#372).

 * Added option `--junc-bed` to load known gene annotations in the BED12
   format. Minimap2 prefers annotated junctions over novel junctions (#197 and
   #348). GTF can be converted to BED12 with `paftools.js gff2bed`.

 * Added option `--sam-hit-only` to suppress unmapped hits in SAM (#377).

 * Added preset `splice:hq` for high-quality CCS or mRNA sequences. It applies
   better scoring and improves the sensitivity to small exons. This preset may
   introduce false small introns, but the overall accuracy should be higher.

This version produces nearly identical alignments to v2.16, except for CIGARs
affected by the bug mentioned above.

(2.17: 5 May 2019, r941)



Release 2.16-r922 (28 February 2019)
------------------------------------

This release is 50% faster for mapping ultra-long nanopore reads at comparable
accuracy. For short-read mapping, long-read overlapping and ordinary long-read
mapping, the performance and accuracy remain similar. This speedup is achieved
with a new heuristic to limit the number of chaining iterations (#324). Users
can disable the heuristic by increasing a new option `--max-chain-iter` to a
huge number.

Other changes to minimap2:

 * Implemented option `--paf-no-hit` to output unmapped query sequences in PAF.
   The strand and reference name columns are both `*` at an unmapped line. The
   hidden option is available in earlier minimap2 but had a different 2-column
   output format instead of PAF.

 * Fixed a bug that leads to wrongly calculated `de` tags when ambiguous bases
   are involved (#309). This bug only affects v2.15.

 * Fixed a bug when parsing command-line option `--splice` (#344). This bug was
   introduced in v2.13.

 * Fixed two division-by-zero cases (#326). They don't affect final alignments
   because the results of the divisions are not used in both case.

 * Added an option `-o` to output alignments to a specified file. It is still
   recommended to use UNIX pipes for on-the-fly conversion or compression.

 * Output a new `rl` tag to give the length of query regions harboring
   repetitive seeds.

Changes to paftool.js:

 * Added a new option to convert the MD tag to the long form of the cs tag.

Changes to mappy:

 * Added the `mappy.Aligner.seq_names` method to return sequence names (#312).

For NA12878 ultra-long reads, this release changes the alignments of <0.1% of
reads in comparison to v2.15. All these reads have highly fragmented alignments
and are likely to be problematic anyway. For shorter or well aligned reads,
this release should produce mostly identical alignments to v2.15.

(2.16: 28 February 2019, r922)



Release 2.15-r905 (10 January 2019)
-----------------------------------

Changes to minimap2:

 * Fixed a rare segmentation fault when option -H is in use (#307). This may
   happen when there are very long homopolymers towards the 5'-end of a read.

 * Fixed wrong CIGARs when option --eqx is used (#266).

 * Fixed a typo in the base encoding table (#264). This should have no
   practical effect.

 * Fixed a typo in the example code (#265).

 * Improved the C++ compatibility by removing "register" (#261). However,
   minimap2 still can't be compiled in the pedantic C++ mode (#306).

 * Output a new "de" tag for gap-compressed sequence divergence.

Changes to paftools.js:

 * Added "asmgene" to evaluate the completeness of an assembly by measuring the
   uniquely mapped single-copy genes. This command learns the idea of BUSCO.

 * Added "vcfpair" to call a phased VCF from phased whole-genome assemblies. An
   earlier version of this script is used to produce the ground truth for the
   syndip benchmark [PMID:30013044].

This release produces identical alignment coordinates and CIGARs in comparison
to v2.14. Users are advised to upgrade due to the several bug fixes.

(2.15: 10 Janurary 2019, r905)



Release 2.14-r883 (5 November 2018)
-----------------------------------

Notable changes:

 * Fixed two minor bugs caused by typos (#254 and #266).

 * Fixed a bug that made minimap2 abort when --eqx was used together with --MD
   or --cs (#257).

 * Added --cap-sw-mem to cap the size of DP matrices (#259). Base alignment may
   take a lot of memory in the splicing mode. This may lead to issues when we
   run minimap2 on a cluster with a hard memory limit. The new option avoids
   unlimited memory usage at the cost of missing a few long introns.

 * Conforming to C99 and C11 when possible (#261).

 * Warn about malformatted FASTA or FASTQ (#252 and #255).

This release occasionally produces base alignments different from v2.13. The
overall alignment accuracy remain similar.

(2.14: 5 November 2018, r883)



Release 2.13-r850 (11 October 2018)
-----------------------------------

Changes to minimap2:

 * Fixed wrongly formatted SAM when -L is in use (#231 and #233).

 * Fixed an integer overflow in rare cases.

 * Added --hard-mask-level to fine control split alignments (#244).

 * Made --MD work with spliced alignment (#139).

 * Replaced musl's getopt with ketopt for portability.

 * Log peak memory usage on exit.

This release should produce alignments identical to v2.12 and v2.11.

(2.13: 11 October 2018, r850)



Release 2.12-r827 (6 August 2018)
---------------------------------

Changes to minimap2:

 * Added option --split-prefix to write proper alignments (correct mapping
   quality and clustered query sequences) given a multi-part index (#141 and
   #189; mostly by @hasindu2008).

 * Fixed a memory leak when option -y is in use.

Changes to mappy:

 * Support the MD/cs tag (#183 and #203).

 * Allow mappy to index a single sequence, to add extra flags and to change the
   scoring system.

Minimap2 should produce alignments identical to v2.11.

(2.12: 6 August 2018, r827)



Release 2.11-r797 (20 June 2018)
--------------------------------

Changes to minimap2:

 * Improved alignment accuracy in low-complexity regions for SV calling. Thank
   @armintoepfer for multiple offline examples.

 * Added option --eqx to encode sequence match/mismatch with the =/X CIGAR
   operators (#156, #157 and #175).

 * When compiled with VC++, minimap2 generated wrong alignments due to a
   comparison between a signed integer and an unsigned integer (#184). Also
   fixed warnings reported by "clang -Wextra".

 * Fixed incorrect anchor filtering due to a missing 64- to 32-bit cast.

 * Fixed incorrect mapping quality for inversions (#148).

 * Fixed incorrect alignment involving ambiguous bases (#155).

 * Fixed incorrect presets: option `-r 2000` is intended to be used with
   ava-ont, not ava-pb. The bug was introduced in 2.10.

 * Fixed a bug when --for-only/--rev-only is used together with --sr or
   --heap-sort=yes (#166).

 * Fixed option -Y that was not working in the previous releases.

 * Added option --lj-min-ratio to fine control the alignment of long gaps
   found by the "long-join" heuristic (#128).

 * Exposed `mm_idx_is_idx`, `mm_idx_load` and `mm_idx_dump` C APIs (#177).
   Also fixed a bug when indexing without reference names (this feature is not
   exposed to the command line).

Changes to mappy:

 * Added `__version__` (#165).

 * Exposed the maximum fragment length parameter to mappy (#174).

Changes to paftools:

 * Don't crash when there is no "cg" tag (#153).

 * Fixed wrong coverage report by "paftools.js call" (#145).

This version may produce slightly different base-level alignment. The overall
alignment statistics should remain similar.

(2.11: 20 June 2018, r797)



Release 2.10-r761 (27 March 2018)
---------------------------------

Changes to minimap2:

 * Optionally output the MD tag for compatibility with existing tools (#63,
   #118 and #137).

 * Use SSE compiler flags more precisely to prevent compiling errors on certain
   machines (#127).

 * Added option --min-occ-floor to set a minimum occurrence threshold. Presets
   intended for assembly-to-reference alignment set this option to 100. This
   option alleviates issues with regions having high copy numbers (#107).

 * Exit with non-zero code on file writing errors (e.g. disk full; #103 and
   #132).

 * Added option -y to copy FASTA/FASTQ comments in query sequences to the
   output (#136).

 * Added the asm20 preset for alignments between genomes at 5-10% sequence
   divergence.

 * Changed the band-width in the ava-ont preset from 500 to 2000. Oxford
   Nanopore reads may contain long deletion sequencing errors that break
   chaining.

Changes to mappy, the Python binding:

 * Fixed a typo in Align.seq() (#126).

Changes to paftools.js, the companion script:

 * Command sam2paf now converts the MD tag to cs.

 * Support VCF output for assembly-to-reference variant calling (#109).

This version should produce identical alignment for read overlapping, RNA-seq
read mapping, and genomic read mapping. We have also added a cook book to show
the variety uses of minimap2 on real datasets. Please see cookbook.md in the
minimap2 source code directory.

(2.10: 27 March 2017, r761)



Release 2.9-r720 (23 February 2018)
-----------------------------------

This release fixed multiple minor bugs.

* Fixed two bugs that lead to incorrect inversion alignment. Also improved the
  sensitivity to small inversions by using double Z-drop cutoff (#112).

* Fixed an issue that may cause the end of a query sequence unmapped (#104).

* Added a mappy API to retrieve sequences from the index (#126) and to reverse
  complement DNA sequences. Fixed a bug where the `best_n` parameter did not
  work (#117).

* Avoided segmentation fault given incorrect FASTQ input (#111).

* Combined all auxiliary javascripts to paftools.js. Fixed several bugs in
  these scripts at the same time.

(2.9: 24 February 2018, r720)



Release 2.8-r672 (1 February 2018)
----------------------------------

Notable changes in this release include:

 * Speed up short-read alignment by ~10%. The overall mapping accuracy stays
   the same, but the output alignments are not always identical to v2.7 due to
   unstable sorting employed during chaining. Long-read alignment is not
   affected by this change as the speedup is short-read specific.

 * Mappy now supports paired-end short-read alignment (#87). Please see
   python/README.rst for details.

 * Added option --for-only and --rev-only to perform alignment against the
   forward or the reverse strand of the reference genome only (#91).

 * Alleviated the issue with undesired diagonal alignment in the self mapping
   mode (#10). Even if the output is not ideal, it should not interfere with
   other alignments. Fully resolving the issue is intricate and may require
   additional heuristic thresholds.

 * Enhanced error checking against incorrect input (#92 and #96).

For long query sequences, minimap2 should output identical alignments to v2.7.

(2.8: 1 February 2018, r672)



Release 2.7-r654 (9 January 2018)
---------------------------------

This release fixed a bug in the splice mode and added a few minor features:

 * Fixed a bug that occasionally takes an intron as a long deletion in the
   splice mode. This was caused by wrong backtracking at the last CIGAR
   operator. The current fix eliminates the error, but it is not optimal in
   that it often produces a wrong junction when the last operator is an intron.
   A future version of minimap2 may improve upon this.

 * Support high-end ARM CPUs that implement the NEON instruction set (#81).
   This enables minimap2 to work on Raspberry Pi 3 and Odroid XU4.

 * Added a C API to construct a minimizer index from a set of C strings (#80).

 * Check scoring specified on the command line (#79). Due to the 8-bit limit,
   excessively large score penalties fail minimap2.

For genomic sequences, minimap2 should give identical alignments to v2.6.

(2.7: 9 January 2018, r654)



Release 2.6-r623 (12 December 2017)
-----------------------------------

This release adds several features and fixes two minor bugs:

 * Optionally build an index without sequences. This helps to reduce the
   peak memory for read overlapping and is automatically applied when
   base-level alignment is not requested.

 * Approximately estimate per-base sequence divergence (i.e. 1-identity)
   without performing base-level alignment, using a MashMap-like method. The
   estimate is written to a new dv:f tag.

 * Reduced the number of tiny terminal exons in RNA-seq alignment. The current
   setting is conservative. Increase --end-seed-pen to drop more such exons.

 * Reduced the peak memory when aligning long query sequences.

 * Fixed a bug that is caused by HPC minimizers longer than 256bp. This should
   have no effect in practice, but it is recommended to rebuild HPC indices if
   possible.

 * Fixed a bug when identifying identical hits (#71). This should only affect
   artifactual reference consisting of near identical sequences.

For genomic sequences, minimap2 should give nearly identical alignments to
v2.5, except the new dv:f tag.

(2.6: 12 December 2017, r623)



Release 2.5-r572 (11 November 2017)
-----------------------------------

This release fixes several bugs and brings a couple of minor improvements:

 * Fixed a severe bug that leads to incorrect mapping coordinates in rare
   corner cases.

 * Fixed underestimated mapping quality for chimeric alignments when the whole
   query sequence contain many repetitive minimizers, and for chimeric
   alignments caused by Z-drop.

 * Fixed two bugs in Python binding: incorrect strand field (#57) and incorrect
   sequence names for Python3 (#55).

 * Improved mapping accuracy for highly overlapping paired ends.

 * Added option -Y to use soft clipping for supplementary alignments (#56).

(2.5: 11 November 2017, r572)



Release 2.4-r555 (6 November 2017)
----------------------------------

As is planned, this release focuses on fine tuning the base algorithm. Notable
changes include

 * Changed the mapping quality scale to match the scale of BWA-MEM. This makes
   minimap2 and BWA-MEM achieve similar sensitivity-specificity balance on real
   short-read data.

 * Improved the accuracy of splice alignment by modeling one additional base
   close to the GT-AG signal. This model is used by default with `-x splice`.
   For SIRV control data, however, it is recommended to add `--splice-flank=no`
   to disable this feature as the SIRV splice signals are slightly different.

 * Tuned the parameters for Nanopore Direct RNA reads. The recommended command
   line is `-axsplice -k14 -uf` (#46).

 * Fixed a segmentation fault when aligning PacBio reads (#47 and #48). This
   bug is very rare but it affects all versions of minimap2. It is also
   recommended to re-index reference genomes created with `map-pb`. For human,
   two minimizers in an old index are wrong.

 * Changed option `-L` in sync with the final decision of hts-specs: a fake
   CIGAR takes the form of `<readLen>S<refLen>N`. Note that `-L` only enables
   future tools to recognize long CIGARs. It is not possible for older tools to
   work with such alignments in BAM (#43 and #51).

 * Fixed a tiny issue whereby minimap2 may waste 8 bytes per candidate
   alignment.

The minimap2 technical note hosted at arXiv has also been updated to reflect
recent changes.

(2.4: 6 November 2017, r555)



Release 2.3-r531 (22 October 2017)
----------------------------------

This release come with many improvements and bug fixes:

 * The **sr** preset now supports paired-end short-read alignment. Minimap2 is
   3-4 times as fast as BWA-MEM, but is slightly less accurate on simulated
   reads.

 * Meticulous improvements to assembly-to-assembly alignment (special thanks to
   Alexey Gurevich from the QUAST team): a) apply a small penalty to matches
   between ambiguous bases; b) reduce missing alignments due to spurious
   overlaps; c) introduce the short form of the `cs` tag, an improvement to the
   SAM MD tag.

 * Make sure gaps are always left-aligned.

 * Recognize `U` bases from Oxford Nanopore Direct RNA-seq (#33).

 * Fixed slightly wrong chaining score. Fixed slightly inaccurate coordinates
   for split alignment.

 * Fixed multiple reported bugs: 1) wrong reference name for inversion
   alignment (#30); 2) redundant SQ lines when multiple query files are
   specified (#39); 3) non-functioning option `-K` (#36).

This release has implemented all the major features I planned five months ago,
with the addition of spliced long-read alignment. The next couple of releases
will focus on fine tuning of the base algorithms.

(2.3: 22 October 2017, r531)



Release 2.2-r409 (17 September 2017)
------------------------------------

This is a feature release. It improves single-end short-read alignment and
comes with Python bindings. Detailed changes include:

 * Added the **sr** preset for single-end short-read alignment. In this mode,
   minimap2 runs faster than BWA-MEM, but is slightly less accurate on
   simulated data sets. Paired-end alignment is not supported as of now.

 * Improved mapping quality estimate with more accurate identification of
   repetitive hits. This mainly helps short-read alignment.

 * Implemented **mappy**, a Python binding for minimap2, which is available
   from PyPI and can be installed with `pip install --user mappy`. Python users
   can perform read alignment without the minimap2 executable.

 * Restructured the indexing APIs and documented key minimap2 APIs in the
   header file minimap.h. Updated example.c with the new APIs. Old APIs still
   work but may become deprecated in future.

This release may output alignments different from the previous version, though
the overall alignment statistics, such as the number of aligned bases and long
gaps, remain close.

(2.2: 17 September 2017, r409)



Release 2.1.1-r341 (6 September 2017)
-------------------------------------

This is a maintenance release that is expected to output identical alignment to
v2.1. Detailed changes include:

 * Support CPU dispatch. By default, minimap2 is compiled with both SSE2 and
   SSE4 based implementation of alignment and automatically chooses the right
   one at runtime. This avoids unexpected errors on older CPUs (#21).

 * Improved Windows support as is requested by Oxford Nanopore (#19). Minimap2
   now avoids variable-length stacked arrays, eliminates alloca(), ships with
   getopt_long() and provides timing functions implemented with Windows APIs.

 * Fixed a potential segmentation fault when specifying -k/-w/-H with
   multi-part index (#23).

 * Fixed two memory leaks in example.c

(2.1.1: 6 September 2017, r341)



Release 2.1-r311 (25 August 2017)
---------------------------------

This release adds spliced alignment for long noisy RNA-seq reads. On a SMRT
Iso-Seq and a Oxford Nanopore data sets, minimap2 appears to outperform
traditional mRNA aligners. For DNA alignment, this release gives almost
identical output to v2.0. Other changes include:

 * Added option `-R` to set the read group header line in SAM.

 * Optionally output the `cs:Z` tag in PAF to encode both the query and the
   reference sequences in the alignment.

 * Fixed an issue where DP alignment uses excessive memory.

The minimap2 technical report has been updated with more details and the
evaluation of spliced alignment:

 * Li, H. (2017). Minimap2: fast pairwise alignment for long nucleotide
   sequences. [arXiv:1708.01492v2](https://arxiv.org/abs/1708.01492v2).

(2.1: 25 August 2017, r311)



Release 2.0-r275 (8 August 2017)
--------------------------------

This release is identical to version 2.0rc1, except the version number. It is
described and evaluated in the following technical report:

 * Li, H. (2017). Minimap2: fast pairwise alignment for long DNA sequences.
   [arXiv:1708.01492v1](https://arxiv.org/abs/1708.01492v1).

(2.0: 8 August 2017, r275)



Release 2.0rc1-r232 (30 July 2017)
----------------------------------

This release improves the accuracy of long-read alignment and added several
minor features.

 * Improved mapping quality estimate for short alignments containing few seed
   hits.

 * Fixed a minor bug that affects the chaining accuracy towards the ends of a
   chain. Changed the gap cost for chaining to reduce false seeding.

 * Skip potentially wrong seeding and apply dynamic programming more frequently.
   This slightly increases run time, but greatly reduces false long gaps.

 * Perform local alignment at Z-drop break point to recover potential inversion
   alignment. Output the SA tag in the SAM format. Added scripts to evaluate
   mapping accuracy for reads simulated with pbsim.

This release completes features intended for v2.0. No major features will be
added to the master branch before the final v2.0.

(2.0rc1: 30 July 2017, r232)



Release r191 (19 July 2017)
---------------------------

This is the first public release of minimap2, an aligner for long reads and
assemblies. This release has a few issues and is generally not recommended for
production uses.

(19 July 2017, r191)
