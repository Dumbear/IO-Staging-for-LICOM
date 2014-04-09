ADIOS_DIR=/home/xuewei/WORK/liangyaxiong/softwares/adios-1.6.0
ADIOS_CONFIG=${ADIOS_DIR}/bin/adios_config

CC = mpicxx
cc = mpicc
LDFLAGS = `${ADIOS_CONFIG} -l` -lskel
CFLAGS  = -g -O2 `${ADIOS_CONFIG} -c`
