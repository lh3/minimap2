The [K8 Javascript shell][k8] is needed to run Javascripts in this directory.
Precompiled k8 binaries for Mac and Linux can be found at the [K8 release
page][k8bin].

* [paf2aln.js](paf2aln.js): convert PAF to [MAF][maf] or BLAST-like output for
  eyeballing. PAF has to be generated with minimap2 option `-S`, which writes
  the aligned sequences to the `cs` tag. An example:
  ```sh
  ../minimap2 -S ../test/MT-*.fa | k8 paf2aln.js /dev/stdin
  ```

* [mapstat.js](mapstat.js): output basic statistics such as the number of
  non-redundant mapped bases, number of split and secondary alignments and
  number of long gaps. This scripts seamlessly works with both SAM and PAF.

* [sim-pbsim.js](sim-pbsim.js): convert reads simulated with [PBSIM][pbsim] to
  FASTA and encode the true mapping positions to read names in a format like
  `S1_33!chr1!225258409!225267761!-`.

* [sim-eval.js](sim-eval.js): evaluate mapping accuracy for FASTA generated
  with [sim-pbsim.js](sim-pbsim.js) or [sim-mason2.js](sim-mason2.js).

* [sam2paf.js](sam2paf.js): convert SAM to PAF.

[k8]: https://github.com/attractivechaos/k8
[k8bin]: https://github.com/attractivechaos/k8/releases
[maf]: https://genome.ucsc.edu/FAQ/FAQformat#format5
[pbsim]: https://github.com/pfaucon/PBSIM-PacBio-Simulator
