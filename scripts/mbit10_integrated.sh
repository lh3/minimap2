#!/bin/bash
make GPU=NV GPUARCH=sm_86 DEBUG=analyze
./minimap2  -t 1 --gpu-chain --gpu-cfg a6000_config.json data/hg38.mmi data/ONT/random_500MBases_90kto100k.fa 