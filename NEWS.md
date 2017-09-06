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
