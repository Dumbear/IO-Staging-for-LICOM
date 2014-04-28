#!/bin/bash

source ../env_config.sh

PROC_WRITER=1
PROC_READER=1
PROC_SERVER=1
let "PROC_ALL=PROC_WRITER+PROC_READER"
DATA_SIZE_X=362
DATA_SIZE_Y=194
DATA_SIZE_Z=30

rm -f conf dataspaces.conf
rm *staged.bp

echo "# Config file for DataSpaces
ndim = 3
dimx = 4096
dimy = 4096
dimz = 4096
max_versions = 2
" > dataspaces.conf

echo "Start DataSpaces server on $PROC_SERVER PEs, -s$PROC_SERVER -c$PROC_ALL"
mpiexec -machinefile server.host -n $PROC_SERVER $DATASPACES -s$PROC_SERVER -c$PROC_ALL >& log.server &

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
mpiexec -machinefile writer.host -n $PROC_WRITER ./licom_skel $DATA_SIZE_X $DATA_SIZE_Y $DATA_SIZE_Z 1 1 1 >& log.licom_skel &

echo "Start reader on $PROC_READER PEs"
mpiexec -machinefile reader.host -n $PROC_READER ./staging 1 3000 $DATA_SIZE_X $DATA_SIZE_Y $DATA_SIZE_Z 1 1 1 >& log.staging &

echo "Wait until all applications exit."
wait
