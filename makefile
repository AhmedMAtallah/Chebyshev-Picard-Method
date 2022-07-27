SRC_DIR := src
IDIR=include 
CC=g++
CFLAGS=-I.

all: mainMCPISHMSerial

SHM.o : $(SRC_DIR)/SHM.c include/SHM.h include/Basic.h include/const.h $(SRC_DIR)/EGM2008.cc
	g++ -c $(SRC_DIR)/SHM.c 

mainMCPISHMSerial : $(SRC_DIR)/mainMCPISHMSerial.cpp $(SRC_DIR)/Basic.cpp $(SRC_DIR)/MCPIIOM.cpp SHM.o 
	g++ $(SRC_DIR)/mainMCPISHMSerial.cpp $(SRC_DIR)/Basic.cpp $(SRC_DIR)/MCPIIOM.cpp SHM.o -o mainMCPISHMSerial -I.

clean:
	rm -f *.o mainMCPISHMSerial



