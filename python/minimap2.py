#!/usr/bin/env python

import sys, getopt
import mappy as mp

def main(argv):
	opts, args = getopt.getopt(argv[1:], "")
	if len(args) < 2:
		print("Usage: minimap2.py <ref.fa>|<ref.mmi> <query.fq>")
		sys.exit(1)
	a = mp.Aligner(args[0]) # load/build index
	if not a: raise Exception("ERROR: failed to load/build index file '{}'".format(args[0]))
	for name, seq, qual in mp.fastx_read(args[1]): # read one sequence
		for h in a.map(seq): # traverse hits
			print('{}\t{}\t{}'.format(name, len(seq), h))

if __name__ == "__main__":
	main(sys.argv)
