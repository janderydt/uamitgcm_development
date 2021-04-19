#!/bin/bash
################################################
# Start a self-resubmitting simulation.
################################################

# ID number for run
JOBNO=05

# clean run directory and link all required files
./prepare_run.sh

# record start times
TIMEQSTART="$(date +%s)"
echo Start-time `date` >> ../run/times

echo DOTSON_$JOBNO
echo $JOBNO
echo $TIMEQSTART
echo $HECACC
# submit the job chain
sbatch --job-name=SMTCPL$JOBNO --account=$HECACC run_repeat_rolling_ckp.slurm
