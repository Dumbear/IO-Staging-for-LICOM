ADIOS_DIR=/home/zouevan/soft/adios-1.6.0
ADIOS_CONFIG=${ADIOS_DIR}/bin/adios_config

CC = mpiicpc
cc = mpiicc
LDFLAGS = `${ADIOS_CONFIG} -l` -lskel
CFLAGS  = -g -O2 `${ADIOS_CONFIG} -c`
