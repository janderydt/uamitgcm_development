#!/bin/sh 
#SBATCH --time=00:20:00
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --tasks-per-node=108
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --reservation=shortqos

# hardwire budget if you wish to over-ride default
#export HECACC=n02-NEI025867
#export HECACC=n02-NEM001660

module load epcc-job-env
export TMPDIR=/work/n02/n02/`whoami`/SCRATCH
export OMP_NUM_THREADS=1

cd $SLURM_SUBMIT_DIR
echo 'MITgcm starts '`date` >> jobs.log
cd ../run



timestep=100
#tottime=124416000
tottime=2177280000
maxntime=181140

rm data

ln -s ../input/* .
ln -s ../build/mitgcmuv .

s=$(grep timeStepNumber pickup.ckptA.meta)
s_sh=$(grep timeStepNumber pickup_shelfice.ckptA.meta)
s_str=$(grep timeStepNumber pickup_streamice.ckptA.meta)
if [ "$s" = "$s_sh" ]; then
	echo "pickups consistent 1"
else
	s=" "
	echo "pickups inconsistent 1"
fi
if [ "$s" = "$s_str" ]; then
	echo "pickups consistent 2"
else
	s=" "
	echo "pickups inconsistent 2"
fi
echo $s
if [ ! -z "$s" ]; then
 IFS=' ' tokens=( $s )
 old=${tokens[3]}
 new=$(echo $old | sed 's/^0*//')
 new2=" niter0=$new"
 new3=" deltaT=$timestep"
 ntime=$(($tottime/$timestep-$new))
 if (( $ntime > $maxntime )); then
  ntime=$maxntime
 fi
 new4=" nTimesteps=$ntime"
 new5=" pickupsuff='ckptA'"
else
 new2=" niter0=0"
 new=0
 new3=" deltaT=$timestep"
 ntime=$(($tottime/$timestep-$new))
 if (( $ntime > $maxntime )); then
  ntime=$maxntime
 fi
 new4=" nTimesteps=$ntime"
 new5="# pickupsuff='ckptA'"
fi
echo $new2
echo $new3
echo $new4
echo $tottime
echo $timestep
echo $new
echo $ntime
sed "s/.*niter0.*/$new2/" data > data.temp
sed "s/.*deltaT.*/$new3/" data.temp > data.temp2
sed "s/.*nTimesteps.*/$new4/" data.temp2 > data.temp3
sed "s/.*pickupsuff.*/$new5/" data.temp3 > data.temp4
mv data.temp4 data

if (( $ntime > 0 )); then
 srun --distribution=block:block --hint=nomultithread ./mitgcmuv > out.txt 2> err.txt
 OUT=$?
else
 OUT=-1
fi

#IFS=' ' tokens=( $s )
#old=${tokens[3]}


#s=$(ls -tr pickup.0*meta | tail -n 1)
s=$(grep timeStepNumber pickup.ckptA.meta)
echo $s
if [ ! -z "$s" ]; then
# IFS='.' tokens=( $s )
# old=${tokens[1]}
 IFS=' ' tokens=( $s )
 old=${tokens[3]}
 new3=$(echo $old | sed 's/^0*//')
else
 new3=0
fi

if (( $new3 > $new )); then
 cp STDOUT.0000 stdout_$new3
 cd ../scripts
# sbatch --job-name=$SLURM_JOB_NAME rput_cirrus_mds.slurm
 if [ $OUT == 0 ]; then 
#  sbatch --job-name=$SLURM_JOB_NAME run_repeat_rolling_ckp.slurm
 fi
else
  echo "no new pickup found"
fi

exit
