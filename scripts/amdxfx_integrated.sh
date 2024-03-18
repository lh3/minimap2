#!/bin/bash
make GPU=AMD DEBUG=analyze
./minimap2  -t 1 --gpu-chain --gpu-cfg gfx1030.json data/hg38.mmi data/ONT/random_500MBases_90kto100k.fa > out.paf
