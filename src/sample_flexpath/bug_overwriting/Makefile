include ../config.mk

all: reader writer

reader: reader.cpp
	$(CXX) $(CFLAGS) -o reader reader.cpp $(LDFLAGS)

writer: writer.cpp
	$(CXX) $(CFLAGS) -o writer writer.cpp $(LDFLAGS)

clean:
	rm -f reader writer
