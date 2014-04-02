ADIOS_DIR=/home/zouevan/soft/adios-1.6.0
ADIOS_CONFIG=${ADIOS_DIR}/bin/adios_config

CC = mpicxx
LDFLAGS = `${ADIOS_CONFIG} -l`
CFLAGS  = -g -O2 `${ADIOS_CONFIG} -c`

all: staging bp_comp

staging: staging.cpp
	$(CC) $(CFLAGS) -o staging staging.cpp $(LDFLAGS)

bp_comp: bp_comp.cpp
	$(CC) $(CFLAGS) -o bp_comp bp_comp.cpp $(LDFLAGS)

clean:
	rm -f staging bp_comp
