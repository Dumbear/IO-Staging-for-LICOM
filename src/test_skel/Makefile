ADIOS_DIR=/home/xuewei/WORK/liangyaxiong/softwares/adios-1.6.0
ADIOS_CONFIG=${ADIOS_DIR}/bin/adios_config

CC = mpicxx
cc = mpicc
LDFLAGS = `${ADIOS_CONFIG} -l` -L${ADIOS_DIR}/lib/skel -lskel
CFLAGS  = -g -O2 `${ADIOS_CONFIG} -c`

all: staging bp_comp licom_skel

staging: staging.cpp
	$(CC) $(CFLAGS) -o staging staging.cpp $(LDFLAGS)

bp_comp: bp_comp.cpp
	$(CC) $(CFLAGS) -o bp_comp bp_comp.cpp $(LDFLAGS)

licom_skel: licom_skel.c
	$(cc) $(CFLAGS) -o licom_skel licom_skel.c $(LDFLAGS)

clean:
	rm -f staging bp_comp
