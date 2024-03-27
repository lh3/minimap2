#!/usr/bin/env python

import sys
import getopt
import mappy as mp

def main(argv):
	opts, args = getopt.getopt(argv[1:], "x:n:m:k:w:r:cM")
	if len(args) < 2:
		print("Usage: minimap2.py [options] <ref.fa>|<ref.mmi> <query.fq>")
		print("Options:")
		print("  -x STR      preset: sr, map-pb, map-ont, asm5, asm10 or splice")
		print("  -n INT      mininum number of minimizers")
		print("  -m INT      mininum chaining score")
		print("  -k INT      k-mer length")
		print("  -w INT      minimizer window length")
		print("  -r INT      band width")
		print("  -c          output the cs tag")
		print("  -M          output the MD tag")
		sys.exit(1)

	preset = min_cnt = min_sc = k = w = bw = None
	out_cs = out_MD = False
	for opt, arg in opts:
		if opt == '-x': preset = arg
		elif opt == '-n': min_cnt = int(arg)
		elif opt == '-m': min_chain_score = int(arg)
		elif opt == '-r': bw = int(arg)
		elif opt == '-k': k = int(arg)
		elif opt == '-w': w = int(arg)
		elif opt == '-c': out_cs = True
		elif opt == '-M': out_MD = True

	a = mp.Aligner(args[0], preset=preset, min_cnt=min_cnt, min_chain_score=min_sc, k=k, w=w, bw=bw)
	if not a: raise Exception("ERROR: failed to load/build index file '{}'".format(args[0]))
	for name, seq, qual in mp.fastx_read(args[1]): # read one sequence
		for h in a.map(seq, cs=out_cs, MD=out_MD): # traverse hits
			print('{}\t{}\t{}'.format(name, len(seq), h))

if __name__ == "__main__":
	main(sys.argv)
