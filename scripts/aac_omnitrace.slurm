#!/bin/bash
#SBATCH --job-name=minimap2         # Job name
#SBATCH --partition=1CN128C8G2H_2IB_MI210_RHEL8               # Specify the MI210 GPU partition or queue
#SBATCH --gres=gpu:1                  # Request 1 GPU
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks-per-node=1           # Number of tasks (processes) per node
#SBATCH --cpus-per-task=16             # Number of CPU cores per task
#SBATCH --mem=500g                      # Memory per node
#SBATCH --time=01:40:00               # Maximum execution time (HH:MM:SS)
#SBATCH --output=slurm_output/sample_sbatch_job.%j.out    # Output file
#SBATCH --error=slurm_output/sample_sbatch_job.%j.err     # Error file

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
# module load rocm/6.0.0
# export AMD_LOG_LEVEL=4
# Replace the following line with the actual command(s) you want to run
cd $MM2_ROOT
make clean
make MICRO_BATCH=5 GPU_CONFIG=aac_config.json SHORT_BLOCK_SIZE=64 LONG_BLOCK_SIZE=1024 MID_BLOCK_SIZE=512 MID_CUT=1 LONG_CUT=100 DEBUG=1 DEBUG_ANALYSIS=1
# omnitrace-sample -PTDH -E all -I rocm-smi -I roctracer -I rocprofiler -I roctx -o omni_output -- ./minimap2 -K 2500000000 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_1GBases_100kto300k.fa
omnitrace-sample -PTDH -E all -I rocm-smi -I roctracer -I rocprofiler -I roctx -o omni_output -- ./minimap2 -K 2000000000 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_5GBases_90kto100k.fa

# ./minimap2 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/reads_4f452f4a-d82a-4580-981b-32d14b997217.fa
# ./minimap2 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_500MBases_200kto300k.fa
# rocprof --stats -o rocprof_output/long_seg.${SLURM_JOB_ID}.csv ./minimap2 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/short_seg_reads_from_1kto10k_distri.fa
# rocprof --stats -o rocprof_output/long_seg.${SLURM_JOB_ID}.csv ./minimap2 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/reads_4f452f4a-d82a-4580-981b-32d14b997217.fa 
# omnitrace-sample -PTDH -E all -I rocm-smi -I roctracer -I rocprofiler -I roctx -o omni_output -- ./minimap2 -K 2000000000 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_500MBases_200kto300k.fa
# omnitrace-sample -PTDH -E all -I rocm-smi -I roctracer -I rocprofiler -I roctx -o omni_output -- ./minimap2 -K 2000000000 -t 1 --max-chain-skip=2147483647 --gpu-chain /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_7GBases_10kto300k.fa
# Optional: You can add post-processing commands here

echo "Exit: $?"

# End of the script
