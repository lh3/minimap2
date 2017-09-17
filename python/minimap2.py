#!/usr/bin/env python

import sys, getopt
import mappy as mp

def readfq(fp): # multi-line fasta/fastq parser
	last = None
	while True:
		if not last:
			for l in fp:
				if l[0] in '>@':
					last = l[:-1]
					break
		if not last: break
		name, seqs, last = last[1:].split()[0], [], None
		for l in fp:
			if l[0] in '@+>':
				last = l[:-1]
				break
			seqs.append(l[:-1])
		if not last or last[0] != '+':
			yield name, ''.join(seqs), None
			if not last: break
		else:
			seq, leng, seqs = ''.join(seqs), 0, []
			for l in fp:
				seqs.append(l[:-1])
				leng += len(l) - 1
				if leng >= len(seq):
					last = None
					yield name, seq, ''.join(seqs);
					break
			if last:
				yield name, seq, None
				break

def main(argv):
	opts, args = getopt.getopt(argv[1:], "")
	if len(args) < 2:
		print("Usage: minimap2.py <ref.fa>|<ref.mmi> <query.fq>")
		sys.exit(1)
	a = mp.Aligner(args[0]) # load/build index
	if not a:
		print("ERROR: failed to load/build index")
		return
	for name, seq, qual in readfq(open(args[1])): # read one sequence
		for h in a.map(seq): # traverse hits
			print('{}\t{}\t{}'.format(name, len(seq), h))

if __name__ == "__main__":
	main(sys.argv)
