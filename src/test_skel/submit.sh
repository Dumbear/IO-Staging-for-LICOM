#!/bin/bash

source ../env_config.sh

PROC_WRITER=1
PROC_READER=1
PROC_SERVER=1
let "PROC_ALL=PROC_WRITER+PROC_READER"
DATA_SIZE_X=64
DATA_SIZE_Y=64
DATA_SIZE_Z=64

rm -f conf dataspaces.conf
rm *staged.bp

echo "# Config file for DataSpaces
ndim = 3
dimx = $DATA_SIZE_X
dimy = $DATA_SIZE_Y
dimz = $DATA_SIZE_Z
max_versions = 2
" > dataspaces.conf

echo "Start DataSpaces server on $PROC_SERVER PEs, -s$PROC_SERVER -c$PROC_ALL"
./bsub_submit.sh 12 $PROC_SERVER 0 $DATASPACES -s$PROC_SERVER -c$PROC_ALL | bsub &> log.server

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
./bsub_submit.sh 12 $PROC_WRITER 0 ./licom_skel $DATA_SIZE_X $DATA_SIZE_Y $DATA_SIZE_Z 1 1 1 | bsub &> log.licom_skel

echo "Start reader on $PROC_READER PEs"
./bsub_submit.sh 12 $PROC_READER 1 ./staging 1 3000 $DATA_SIZE_X $DATA_SIZE_Y $DATA_SIZE_Z 1 1 1 | bsub &>  log.staging

echo "Wait until all applications exit."
wait
