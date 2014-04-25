ADIOS_DIR=/home/zouevan/soft/adios-1.6.0
ADIOS_CONFIG=${ADIOS_DIR}/bin/adios_config

CC = mpiicpc -cxx="/home/software/intel_updata/composer_xe_2013_sp1.2.144/bin/intel64/icpc"
cc = mpiicc -cc="/home/software/intel_updata/composer_xe_2013_sp1.2.144/bin/intel64/icc"
LDFLAGS = `${ADIOS_CONFIG} -l` -lskel
CFLAGS  = -g -O2 `${ADIOS_CONFIG} -c`
