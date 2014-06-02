#!/bin/bash

source ../env_config.sh

PROC_WRITER=4
PROC_READER=1

rm *_info.txt *_ready.txt

echo "Start writer on $PROC_WRITER PEs"
mpiexec -host cn001 -n $PROC_WRITER ./writer >& log.writer &

echo "Start reader on $PROC_READER PEs"
mpiexec -host cn001 -n $PROC_READER ./reader >& log.reader &

echo "Wait until all applications exit."
wait
