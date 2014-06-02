#!/bin/bash

source ../env_config.sh

PROC_WRITER=4
PROC_READER=1

rm *_info.txt *_ready.txt

echo "Start writer on $PROC_WRITER PEs"
mpiexec -host cn001 -n $PROC_WRITER ./writer1 >& log.writer1 &

echo "Start reader on $PROC_READER PEs"
mpiexec -host cn001 -n $PROC_READER ./reader

echo "Wait until all applications exit."
wait

rm *_info.txt *_ready.txt

echo "Start writer on $PROC_WRITER PEs"
mpiexec -host cn001 -n $PROC_WRITER ./writer2 >& log.writer2 &

echo "Start reader on $PROC_READER PEs"
mpiexec -host cn001 -n $PROC_READER ./reader

echo "Wait until all applications exit."
wait
