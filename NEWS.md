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
