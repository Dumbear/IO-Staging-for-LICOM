ADIOS_DIR=/home/liangyaxiong/softwares/ADIOS-master
ADIOS_CONFIG=${ADIOS_DIR}/bin/adios_config

CXX = mpiicpc -cxx="/home/software/intel_updata/composer_xe_2013_sp1.2.144/bin/intel64/icpc"
CC = mpiicc -cc="/home/software/intel_updata/composer_xe_2013_sp1.2.144/bin/intel64/icc"
LDFLAGS = `${ADIOS_CONFIG} -l` -lskel
CFLAGS  = -g -O2 `${ADIOS_CONFIG} -c`
