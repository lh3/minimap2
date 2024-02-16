ncu --export "ncu/report" --force-overwrite --target-processes all --kernel-name-base function --kernel-name score_generation_long --launch-skip-before-match 0 --section ComputeWorkloadAnalysis --section InstructionStats --section LaunchStats --section MemoryWorkloadAnalysis --section MemoryWorkloadAnalysis_Chart --section MemoryWorkloadAnalysis_Tables --section Occupancy --section SchedulerStats --section SourceCounters --section SpeedOfLight --section SpeedOfLight_RooflineChart --section WarpStateStats --sampling-interval auto --sampling-max-passes 5 --sampling-buffer-size 33554432 --profile-from-start 1 --cache-control all --clock-control base --apply-rules yes --import-source yes --source-folders gpu --check-exit-code yes --replay-mode application ./minimap2 data/hg38.mmi data/random_100MBases_1kto300k.fa -t 4 --max-chain-skip=2147483647 --gpu-chain