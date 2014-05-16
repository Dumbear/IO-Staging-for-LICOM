#!/bin/bash

source ../env_config.sh

PROC_WRITER=1
PROC_READER=1
PROC_SERVER=1
let "PROC_ALL=PROC_WRITER+PROC_READER"
DATA_SIZE_X=362
DATA_SIZE_Y=194
DATA_SIZE_Z=30

rm *_info.txt *_ready.txt
rm *staged.bp

echo "Start writer on $PROC_WRITER PEs"
mpiexec -machinefile writer.host -n $PROC_WRITER ./licom_skel $DATA_SIZE_X $DATA_SIZE_Y $DATA_SIZE_Z 1 1 1 >& log.licom_skel &

echo "Start reader on $PROC_READER PEs"
mpiexec -machinefile reader.host -n $PROC_READER ./staging 1 3000 $DATA_SIZE_X $DATA_SIZE_Y $DATA_SIZE_Z 1 1 1 >& log.staging &

echo "Wait until all applications exit."
wait
