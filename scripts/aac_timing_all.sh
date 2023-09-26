#!/bin/bash
WORKSPACE_DIR=/shared/prod/home/liuxs/bioinfo/minimap2
EXE_PATH=$WORKSPACE_DIR
CONFIG_PATH=$WORKSPACE_DIR
# DATA_PATH=/shared/prod/home/liuxs/bioinfo/Profile_mm2/data
DATA_PATH=/shareddata/umich_folder/data/ONT
# NOTE: currently `-t n_thread`, the n_thread must be equal to the num_streams in gpu_config.json
N_THREAD=$(sed -n 's/.*"num_streams": \([0-9]*\).*/\1/p' ${CONFIG_PATH}/gpu_config.json)
# Array of LONG_BLOCK_SIZE values
MID_BLOCK_SIZES=( 512 )
MID_CUTS=( 1 )
LONG_CUTS=( 100 )
GPU_CONFIGS=( gpu_config.json )
# DATA_SETS=( 50kto100k )
# DATA_SETS=( 1kto5k 1kto10k 1kto50k 1kto300k 200kto300k 1kto20k 1kto30k 1kto70k 1kto200k 10kto50k 10kto100k 50kto100k)
# DATA_SETS=( 1kto300k 50kto300k 100kto300k 150kto300k 200kto300k 250kto300k 1kto200k 20kto200k 50kto200k 70kto200k 100kto200k 130kto200k 150kto200k 170kto200k)
DATA_SETS=( 1kto5k 9kto10k 10kto20k 20kto30k 40kto50k 90kto100k 110kto120k 140kto150k 180kto200k 200kto250k 200kto300k )

# Iterate over LONG_BLOCK_SIZES array
for MID_BLOCK_SIZE in "${MID_BLOCK_SIZES[@]}"
do
    for MID_CUT in "${MID_CUTS[@]}"
    do
        for LONG_CUT in "${LONG_CUTS[@]}"
        do
            echo "Executing with MID_BLOCK_SIZE=${MID_BLOCK_SIZE} MID_CUT=${MID_CUT} LONG_CUT=${LONG_CUT}"
            make clean 
            make GPU_CONFIG=${GPU_CONFIG} SHORT_BLOCK_SIZE=64 MID_BLOCK_SIZE=${MID_BLOCK_SIZE} LONG_BLOCK_SIZE=1024 MID_CUT=${MID_CUT} LONG_CUT=${LONG_CUT} 
              
            for GPU_CONFIG in "${GPU_CONFIGS[@]}"
            do
                for DATA_SET in "${DATA_SETS[@]}"
                do
                    QUERY_FILE=$DATA_PATH/random_500MBases_${DATA_SET}.fa
                    echo "Executing on dataset ${DATA_SET}"    
                    filename="profile_output/data-${DATA_SET}_profile_${N_THREAD}_midblk-${MID_BLOCK_SIZE}_cut-${LONG_CUT}"                
                    ${EXE_PATH}/minimap2 ${DATA_PATH}/hg38.mmi -t ${N_THREAD} --max-chain-skip=2147483647 --gpu-chain ${QUERY_FILE} > test.out 2> $filename
                done
            done
        done
    done
done