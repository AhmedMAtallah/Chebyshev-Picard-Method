all: mainMCPISHMSerial

SHM.o : SHM.c SHM.h Basic.h const.h EGM2008.cc
	g++ -c SHM.c 

mainMCPISHMSerial : mainMCPISHMSerial.cpp Basic.cpp MCPIIOM.cpp SHM.o 
	g++ mainMCPISHMSerial.cpp Basic.cpp MCPIIOM.cpp SHM.o -o mainMCPISHMSerial -I.

clean:
	rm -f *.o mainMCPISHMSerial



