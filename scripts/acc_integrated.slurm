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

# random data test
# export AMD_LOG_LEVEL=4
make clean
# debug level: N/A info analyze verbose
make GPU=AMD DEBUG=analyze
# ./minimap2 -K 2000000000 -t 1 --gpu-chain --gpu-cfg /shared/prod/home/liuxs/bioinfo/minimap2/mi210_over50k.json /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_5GBases_90kto100k.fa > out.paf
# echo "Exit: $?"
./minimap2 -K 2000000000 -t 1 --gpu-chain --gpu-cfg /shared/prod/home/liuxs/bioinfo/minimap2/mi210_below50k.json /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_5GBases_20kto30k.fa > out.paf
echo "Exit: $?"
./minimap2 -K 2000000000 -t 1 --gpu-chain --gpu-cfg /shared/prod/home/liuxs/bioinfo/minimap2/mi210_below50k.json /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_5GBases_60kto70k.fa > out.paf
echo "Exit: $?"
./minimap2 -K 2000000000 -t 1 --gpu-chain --gpu-cfg /shared/prod/home/liuxs/bioinfo/minimap2/mi210_over50k.json /shareddata/umich_folder/data/ONT/hg38.mmi /shareddata/umich_folder/data/ONT/random_5GBases_60kto70k.fa > out.paf
echo "Exit: $?"
# Optional: You can add post-processing commands here

# End of the script
