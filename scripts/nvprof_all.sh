#!/bin/bash

# Array of LONG_BLOCK_SIZE values
MID_BLOCK_SIZES=( 128 256 512 768 1024 )
MID_CUTS=( 1 )
LONG_CUTS=( 50 60 70 )

# Iterate over LONG_BLOCK_SIZES array
for MID_BLOCK_SIZE in "${MID_BLOCK_SIZES[@]}"
do
    for MID_CUT in "${MID_CUTS[@]}"
    do
        for LONG_CUT in "${LONG_CUTS[@]}"
        do
            echo "Executing with MID_BLOCK_SIZE=${MID_BLOCK_SIZE} MID_CUT=${MID_CUT} LONG_CUT=${LONG_CUT}"

            # Clean the project
            make clean

            # Build with specific configurations
            make GPU=NVCC GPU_CONFIG=gpu_config.json SHORT_BLOCK_SIZE=64 MID_BLOCK_SIZE=${MID_BLOCK_SIZE} LONG_BLOCK_SIZE=1024 MID_CUT=${MID_CUT} LONG_CUT=${LONG_CUT} 

            # # Run nsys nvprof with the given command # data/random_500MBases_200kto300k.fa
            # read by `nsys stats --report cuda_gpu_kern_sum nsys/nvprofxxx_xxx_xxx.nsys-rep` 
            nsys nvprof -o nsys/nvprof_mid${MID_BLOCK_SIZE}_${MID_CUT}_${LONG_CUT} ./minimap2 data/hg38.mmi data/random_500MBases_200kto300k.fa -t 1 --max-chain-skip=2147483647 --gpu-chain 2> nsys/runtime_report_mid${MID_BLOCK_SIZE}_${MID_CUT}_${LONG_CUT}.txt
            # nsys stats --format csv --report cuda_gpu_kern_sum nsys/nvprofxxx_xxx_xxx.nsys-rep
            echo "Done with MID_BLOCK_SIZE=${MID_BLOCK_SIZE} MID_CUT=${MID_CUT} LONG_CUT=${LONG_CUT}"
        done
    done
done
