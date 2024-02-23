#!/bin/bash

make clean
make MICRO_BATCH=4 GPU=AMD GPU_CONFIG=gtx1030.json SHORT_BLOCK_SIZE=32 LONG_BLOCK_SIZE=1024 MID_BLOCK_SIZE=512 MID_CUT=1 LONG_CUT=100 DEBUG_ANALYSIS=1
rocgdb --args ./minimap2  -t 1 --gpu-chain --max-chain-skip=2147483647 data/hg38.mmi data/ONT/random_500MBases_90kto100k.fa 