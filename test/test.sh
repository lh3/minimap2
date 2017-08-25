#!/bin/bash
# test build script and a few simple test cases 

# to die on any error (non-zero return)
set -e 

# test build
echo "0. Test build"
make && echo "passed"

# test cases removing files generated if test passes
echo "1. Testing simple mapping"
./minimap2 -ax map10k test/MT-human.fa test/MT-orang.fa > test_simple.sam && echo "passed" 
rm test_simple.sam

echo "2. Testing creating an index first and then mapping"
./minimap2 -x map10k -d test_MT-human.mmi test/MT-human.fa
./minimap2 -ax map10k test_MT-human.mmi test/MT-orang.fa > test_index.sam && echo "passed"
rm test_MT-human.mmi test_index.sam
