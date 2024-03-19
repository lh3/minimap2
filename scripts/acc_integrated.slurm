#!/bin/bash
#SBATCH --job-name=minimap2         # Job name
#SBATCH --partition=1CN128C8G2H_2IB_MI210_RHEL8               # Specify the MI210 GPU partition or queue
#SBATCH --gres=gpu:1                  # Request 1 GPU
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks-per-node=1           # Number of tasks (processes) per node
#SBATCH --cpus-per-task=16             # Number of CPU cores per task
#SBATCH --mem=500g                      # Memory per node
#SBATCH --time=02:40:00               # Maximum execution time (HH:MM:SS)
#SBATCH --output=slurm_output/sample_sbatch_job.%j.out    # Output file
#SBATCH --error=slurm_output/sample_sbatch_job.%j.err     # Error file

#salloc --partition=1CN128C8G2H_2IB_MI210_RHEL8 --gres=gpu:1 --nodes=1 --cpus-per-task=8 --time=01:00:00

export MM2_ROOT=$HOME/minimap2

# Load necessary modules (if required)
source /etc/profile.d/modules.sh
scl enable gcc-toolset-11 bash
module unuse /shared/apps/modules/ubuntu/modulefiles
module use /shared/apps/modules/rhel8/modulefiles
module unuse /shared/apps/modules/rhel9/modulefiles
module unuse /shared/apps/modules/sles15sp4/modulefiles
module unuse /shared/apps/modules/centos8/modulefiles
module unuse /shared/apps/modules/rocky9/modulefiles

module purge
module load rocm-5.7.1
# module load rocm-6.0.0
# Replace the following line with the actual command(s) you want to run
cd $MM2_ROOT
# SOL exp 
# make clean
# make MICRO_BATCH=2 GPU_CONFIG=aac_config.json SHORT_BLOCK_SIZE=64 LONG_BLOCK_SIZE=1024 MID_BLOCK_SIZE=512 MID_CUT=1 LONG_CUT=100 DEBUG=1 DEBUG_ANALYSIS=1
# ./minimap2 -K 2000000000 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/reads_4f452f4a-d82a-4580-981b-32d14b997217.fa
# long_read_600M.fa 

# random data test
# export AMD_LOG_LEVEL=4
# make clean
# make MICRO_BATCH=5 GPU_CONFIG=aac_config.json SHORT_BLOCK_SIZE=64 LONG_BLOCK_SIZE=1024 MID_BLOCK_SIZE=512 MID_CUT=1 LONG_CUT=40 DEBUG=1 DEBUG_ANALYSIS=1

# ./minimap2 -K 2000000000 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_500MBases_1kto5k.fa
# echo "Exit: $?"
# ./minimap2 -K 2000000000 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_500MBases_10kto20k.fa
# echo "Exit: $?"
# ./minimap2 -K 2000000000 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_500MBases_40kto50k.fa
# echo "Exit: $?"
# ./minimap2 -K 2000000000 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_5GBases_10kto300k.fa
# echo "Exit: $?"
# ./minimap2 -K 2000000000 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_500MBases_90kto100k.fa
# echo "Exit: $?"
# ./minimap2 -K 2000000000 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_500MBases_140kto150k.fa
# echo "Exit: $?"
# ./minimap2 -K 2000000000 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_4GBases_10kto300k.fa 


# ./minimap2 -K 2000000000 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/long_read_3G.fa 
# ./minimap2 -K 2000000000 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/long_read_600M.fa 
# ./minimap2 -K 2000000000 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_5GBases_90kto100k.fa
# ./minimap2 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_500MBases_90kto100k.fa
# rocgdb --args ./minimap2 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi random_250reads_90kto100k.fa
# ./minimap2 -K 2000000000 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_500MBases_1kto300k.fa
# ./minimap2 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/short_seg_reads_from_1kto10k_distri.fa

# Optional: You can add post-processing commands here


# Iterate from 2 to 10 for Mid seg length
for ((i = 20; i <= 20; i+=30))
do
    make clean
    make MICRO_BATCH=2 GPU_CONFIG=aac_config.json SHORT_BLOCK_SIZE=64 LONG_BLOCK_SIZE=1024 MID_BLOCK_SIZE=512 MID_CUT=3 LONG_CUT=$i DEBUG=1 DEBUG_ANALYSIS=1
    ./minimap2 -K 2000000000 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_500MBases_40kto50k.fa
    echo "Exit: $?"
    mv range_dis.csv data/range_dis_l$i.csv
    mv mid_range_dis.csv data/mid_range_dis_l$i.csv
    mv long_range_dis.csv data/long_range_dis_l$i.csv
done

# End of the script
