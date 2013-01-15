CC=g++
CFLAGS= -std=c++0x -c -Wall

all: gil

gil: hazard.o pnet.o main.o
	$(CC) main.o pnet.o hazard.o -o gil

hazard.o: hazard.cc hazard.hh
	$(CC) $(CFLAGS) hazard.cc

pnet.o: pnet.cc pnet.hh hazard.hh
	$(CC) $(CFLAGS) pnet.cc hazard.o	

main.o: main.cc
	$(CC) $(CFLAGS) main.cc pnet.o hazard.o

clean:
	rm -rf *o gil