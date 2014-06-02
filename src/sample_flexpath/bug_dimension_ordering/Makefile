include ../config.mk

all: reader writer1 writer2

reader: reader.cpp
	$(CXX) $(CFLAGS) -o reader reader.cpp $(LDFLAGS)

writer1: writer1.cpp
	$(CXX) $(CFLAGS) -o writer1 writer1.cpp $(LDFLAGS)

writer2: writer2.cpp
	$(CXX) $(CFLAGS) -o writer2 writer2.cpp $(LDFLAGS)

clean:
	rm -f reader writer1 writer2
