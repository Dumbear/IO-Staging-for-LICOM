#!/bin/bash

source ../env_config.sh

PROC_WRITER=1
PROC_SERVER=1
WRITE_SIZE=1
let "TOTAL_SIZE=WRITE_SIZE*1024*1024"

rm -f conf dataspaces.conf

echo "# Config file for DataSpaces
ndim = 3
dimx = $TOTAL_SIZE
dimy = 1
dimz = 1
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
bsub -q normal -a intelmpi -n $PROC_WRITER -o log.writer mpirun.lsf ./writer $WRITE_SIZE

echo "Wait until all applications exit."
wait
