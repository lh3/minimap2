#!/bin/bash

# Array of LONG_BLOCK_SIZE values
MID_BLOCK_SIZES=( 512 256 )
MID_CUTS=( 1 )
LONG_CUTS=( 50 100 )
GPU_CONFIGS=( gpu_config1.json )
DATA_SETS=( 1kto300k )

# Iterate over LONG_BLOCK_SIZES array
for DATA_SET in "${DATA_SETS[@]}"
do
    for MID_BLOCK_SIZE in "${MID_BLOCK_SIZES[@]}"
    do
        for MID_CUT in "${MID_CUTS[@]}"
        do
            for LONG_CUT in "${LONG_CUTS[@]}"
            do
                for GPU_CONFIG in "${GPU_CONFIGS[@]}"
                do
                    echo "Executing with MID_BLOCK_SIZE=${MID_BLOCK_SIZE} MID_CUT=${MID_CUT} LONG_CUT=${LONG_CUT}"

                    # Clean the project
                    make GPU=NVCC clean

                    # Build with specific configurations
                    make GPU=NVCC GPU_CONFIG=${GPU_CONFIG} SHORT_BLOCK_SIZE=64 MID_BLOCK_SIZE=${MID_BLOCK_SIZE} LONG_BLOCK_SIZE=512 MID_CUT=${MID_CUT} LONG_CUT=${LONG_CUT} 
                    
                    # # Run nsys ncu with the given command # data/random_500MBases_200kto300k.fa
                    # read by `ncu --import ncu_midXXX_X_XXX_gpu_configX.json_XXXktoXXXk_report.ncu-rep --section SpeedOfLight --kernel-name score_generation_mid -c 1 --csv > test.csv` 
                    ncu --export "ncu/ncu32_mid${MID_BLOCK_SIZE}_${MID_CUT}_${LONG_CUT}_${GPU_CONFIG}_${DATA_SET}_report" --force-overwrite --target-processes all --kernel-name-base function --kernel-name regex:score_generation_mid --launch-count 1 --launch-skip-before-match 0 --section ComputeWorkloadAnalysis --section InstructionStats --section LaunchStats --section MemoryWorkloadAnalysis --section MemoryWorkloadAnalysis_Chart --section MemoryWorkloadAnalysis_Tables --section Occupancy --section SchedulerStats --section SourceCounters --section SpeedOfLight --section SpeedOfLight_RooflineChart --section WarpStateStats --sampling-interval auto --sampling-max-passes 5 --sampling-buffer-size 33554432 --profile-from-start 1 --cache-control all --clock-control base --apply-rules yes --import-source yes --source-folders gpu --check-exit-code yes ./minimap2 data/hg38.mmi data/random_500MBases_${DATA_SET}.fa -t 1 --max-chain-skip=2147483647 --gpu-chain # 2> ncu/runtime_report_mid${MID_BLOCK_SIZE}_${MID_CUT}_${LONG_CUT}_${DATA_SET}.txt
                    echo "Done with MID_BLOCK_SIZE=${MID_BLOCK_SIZE} MID_CUT=${MID_CUT} LONG_CUT=${LONG_CUT}"
                done
            done
        done
    done
done