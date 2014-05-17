include ../config.mk

all: staging licom_skel

staging: staging.cpp
	$(CXX) $(CFLAGS) -o staging staging.cpp $(LDFLAGS)

licom_skel: licom_skel.c
	$(CC) $(CFLAGS) -o licom_skel licom_skel.c $(LDFLAGS)

clean:
	rm -f staging licom_skel
