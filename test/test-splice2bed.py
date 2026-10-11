#!/usr/bin/env python3
"""Run with K8=/path/to/k8 python3 test/test-splice2bed.py.

PAFTOOLS may point to another paftools.js to check an unpatched version.
"""
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest


K8 = os.environ.get("K8", "k8")
PAFTOOLS = os.environ.get(
    "PAFTOOLS", str(Path(__file__).resolve().parents[1] / "misc" / "paftools.js")
)


def alignment(cigar, fmt="sam", strand="+", name="read1", secondary=False):
    operations = [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]
    reference_length = sum(n for n, op in operations if op in "MDN=X")
    query_length = sum(n for n, op in operations if op in "MIS=X")
    if fmt == "sam":
        flag = (16 if strand == "-" else 0) | (256 if secondary else 0)
        fields = [name, flag, "chr1", 101, 60, cigar, "*", 0, 0, "*", "*"]
    else:
        fields = [name, query_length, 0, query_length, strand, "chr1", 10000,
                  100, 100 + reference_length, 30, 30, 60, "cg:Z:" + cigar]
        if not secondary:
            fields.append("s2:i:0")
    return "\t".join(map(str, fields)) + "\n"


def convert(data, options=(), from_file=False):
    command = [K8, PAFTOOLS, "splice2bed", *options]
    if from_file:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "alignments.sam"
            path.write_text(data)
            result = subprocess.run(command + [str(path)], text=True,
                                    capture_output=True, check=True)
    else:
        result = subprocess.run(command + ["-"], input=data, text=True,
                                capture_output=True, check=True)
    return [line.split("\t") for line in result.stdout.splitlines()]


class Splice2BedTests(unittest.TestCase):
    def assert_blocks(self, record, end, sizes, starts, strand="+"):
        self.assertEqual(record, ["chr1", "100", str(end), "read1", "1000", strand,
                                  "100", str(end), "0,128,255", str(len(sizes)),
                                  ",".join(map(str, sizes)) + ",",
                                  ",".join(map(str, starts)) + ","])

    def test_equivalent_match_operations(self):
        for fmt in ("sam", "paf"):
            for strand in ("+", "-"):
                for cigar in ("10M100N20M", "10=100N20=", "8=2X100N19=1X"):
                    with self.subTest(fmt=fmt, strand=strand, cigar=cigar):
                        records = convert(alignment(cigar, fmt, strand))
                        self.assertEqual(len(records), 1)
                        self.assert_blocks(records[0], 230, [10, 20], [0, 110], strand)

    def test_adjacent_introns_have_no_empty_internal_block(self):
        for fmt in ("sam", "paf"):
            for cigar in ("10M100N3I100N20M", "10=100N3I100N20=",
                          "10M100N100N20M", "10M100N3I100N2I50N20M"):
                with self.subTest(fmt=fmt, cigar=cigar):
                    extra = 50 if "50N" in cigar else 0
                    records = convert(alignment(cigar, fmt))
                    self.assert_blocks(records[0], 330 + extra, [10, 20], [0, 210 + extra])

    def test_insertions_clipping_padding_and_deletions(self):
        for fmt in ("sam", "paf"):
            with self.subTest(fmt=fmt):
                cigar = "2H3S8=2X2I3D100N5M1P15=4S1H"
                records = convert(alignment(cigar, fmt))
                self.assert_blocks(records[0], 233, [13, 20], [0, 113])

    def test_unspliced_records(self):
        for fmt in ("sam", "paf"):
            for cigar in ("30M", "30=", "28=2X", "10M2D18M"):
                with self.subTest(fmt=fmt, cigar=cigar):
                    records = convert(alignment(cigar, fmt))
                    self.assert_blocks(records[0], 130, [30], [0])

    def test_headers_unmapped_and_separate_read_groups(self):
        data = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10000\n"
        data += "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n"
        data += alignment("10=100N20=") + alignment("30M", name="read2")
        records = convert(data, from_file=True)
        self.assertEqual([record[3] for record in records], ["read1", "read2"])
        self.assert_blocks(records[0], 230, [10, 20], [0, 110])

    def test_secondary_filtering(self):
        for fmt in ("sam", "paf"):
            with self.subTest(fmt=fmt):
                data = alignment("10=100N20=", fmt)
                data += alignment("30=", fmt, secondary=True)
                self.assertEqual(len(convert(data)), 1)
                records = convert(data, options=("-m",))
                self.assertEqual(len(records), 2)
                self.assertEqual(records[1][8], "0,192,0")
                self.assertEqual(records[1][10:], ["30,", "0,"])


if __name__ == "__main__":
    unittest.main()
