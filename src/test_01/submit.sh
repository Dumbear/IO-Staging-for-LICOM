#!/bin/bash

source ../env_config.sh

PROC_WRITER=1
PROC_SERVER=1
DATA_SIZE_X=1024
DATA_SIZE_Y=1024
DATA_SIZE_Z=1024

rm -f conf dataspaces.conf

echo "# Config file for DataSpaces
ndim = 3
dimx = $DATA_SIZE_X
dimy = $DATA_SIZE_Y
dimz = $DATA_SIZE_Z
max_versions = 2
" > dataspaces.conf

echo "Start DataSpaces server on $PROC_SERVER PEs, -s$PROC_SERVER -c$PROC_WRITER"
bsub -q normal -a intelmpi -n $PROC_SERVER -o log.server mpirun.lsf $DATASPACES -s$PROC_SERVER -c$PROC_WRITER

sleep 1s
while [ ! -f conf ]; do
    echo "File conf is not yet available from server, sleep more..."
    sleep 1s
done
sleep 3s

while read line; do
    export set "$line"
done < conf

echo "DataSpaces IDs: P2TNID = $P2TNID, P2TPID = $P2TPID"

echo "Start writer on $PROC_WRITER PEs"
bsub -q normal -a intelmpi -n $PROC_WRITER -o log.writer mpirun.lsf ./writer $DATA_SIZE_X $DATA_SIZE_Y $DATA_SIZE_Z 1 1 1

echo "Wait until all applications exit."
wait
