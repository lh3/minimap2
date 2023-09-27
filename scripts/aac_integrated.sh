#!/bin/bash
WORKSPACE_DIR=/shared/prod/home/liuxs/bioinfo/minimap2
EXE_PATH=$WORKSPACE_DIR
CONFIG_PATH=$WORKSPACE_DIR
DATA_PATH=/shared/prod/home/liuxs/bioinfo/Profile_mm2/data
# NOTE: currently `-t n_thread`, the n_thread must be equal to the num_streams in gpu_config.json
N_THREAD=$(sed -n 's/.*"num_streams": \([0-9]*\).*/\1/p' ${CONFIG_PATH}/gpu_config.json)
# Array of LONG_BLOCK_SIZE values
MID_BLOCK_SIZES=( 64 )
MID_CUTS=( 1 )
LONG_CUTS=( 50 100 )
GPU_CONFIGS=( gpu_config.json )
DATA_SETS=( 1kto300k 200kto300k )

# Iterate over LONG_BLOCK_SIZES array
for DATA_SET in "${DATA_SETS[@]}"
do
    QUERY_FILE=$DATA_PATH/random_500MBases_${DATA_SET}.fa
    for MID_BLOCK_SIZE in "${MID_BLOCK_SIZES[@]}"
    do
        for MID_CUT in "${MID_CUTS[@]}"
        do
            for LONG_CUT in "${LONG_CUTS[@]}"
            do
                for GPU_CONFIG in "${GPU_CONFIGS[@]}"
                do
                    echo "Executing with MID_BLOCK_SIZE=${MID_BLOCK_SIZE} MID_CUT=${MID_CUT} LONG_CUT=${LONG_CUT}"

                    make clean 

                    make GPU_CONFIG=${GPU_CONFIG} SHORT_BLOCK_SIZE=64 MID_BLOCK_SIZE=${MID_BLOCK_SIZE} LONG_BLOCK_SIZE=1024 MID_CUT=${MID_CUT} LONG_CUT=${LONG_CUT} 
                                        
                    omniperf profile -n integrated_${MID_BLOCK_SIZE}_${LONG_CUT}_${DATA_SET}_report --device 0 -- ${EXE_PATH}/minimap2 ${DATA_PATH}/hg38.mmi -t ${N_THREAD} --max-chain-skip=2147483647 --gpu-chain ${QUERY_FILE} > test_logs.out
                done
            done
        done
    done
done