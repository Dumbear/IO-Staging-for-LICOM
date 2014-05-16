include ../config.mk

all: staging licom_skel

staging: staging.cpp
	$(CC) $(CFLAGS) -o staging staging.cpp $(LDFLAGS)

licom_skel: licom_skel.c
	$(cc) $(CFLAGS) -o licom_skel licom_skel.c $(LDFLAGS)

clean:
	rm -f staging licom_skel
