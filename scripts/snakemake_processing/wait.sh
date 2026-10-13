#!/bin/bash

jobid=$(sbatch --wait --parsable "$@")
echo "Submitted job $jobid"
exit $?


# while squeue -j "$jobid" 2>/dev/null | grep -q "$jobid"; do
#     sleep 30
# done