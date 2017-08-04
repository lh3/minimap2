The [K8 Javascript engine][k8] is needed to run Javascripts in this directory.
Precompiled k8 binaries for Mac and Linux can be found at the [K8 release
page][k8].

* [paf2aln.js](paf2aln.js): convert PAF to [MAF][maf] or BLAST-like output for
  eyeballing. PAF has to be generated with minimap2 option `-S`.

* [mapstat.js](mapstat.js): output basic statistics such as the number of
  non-redundant mapped bases, number of split and secondary alignments and
  number of long gaps. This scripts seamlessly works with both SAM and PAF.

* [pbsim2fa.js](pbsim2fa.js): convert reads simulated with [PBSIM][pbsim] to
  FASTA and encode the true mapping positions to read names.

* [pbsim-eval.js](pbsim-eval.js): evaluate mapping accuracy for FASTA generated
  with [pbsim2fa.js](pbsim2fa.js). This script only works with PAF. For SAM,
  please run the following instead:
  ```sh
  k8 sam2paf.js -p aln.sam | k8 pbsim-eval.js /dev/stdin
  ```

* [sam2paf.js](sam2paf.js): convert SAM to PAF.

[k8]: https://github.com/attractivechaos/k8
[k8bin]: https://github.com/attractivechaos/k8/releases
[maf]: https://genome.ucsc.edu/FAQ/FAQformat#format5
[pbsim]: https://github.com/pfaucon/PBSIM-PacBio-Simulator
