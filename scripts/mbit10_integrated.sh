#!/bin/bash

make clean
make MICRO_BATCH=4 GPU=NV GPU_CONFIG=a6000_config.json SHORT_BLOCK_SIZE=32 LONG_BLOCK_SIZE=1024 MID_BLOCK_SIZE=512 MID_CUT=1 LONG_CUT=100 DEBUG_ANALYSIS=1
./minimap2  -t 1 --max-chain-skip=2147483647 --gpu-chain data/hg38.mmi data/random_500MBases_90kto100k.fa 