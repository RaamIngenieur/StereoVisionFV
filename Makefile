CC=g++
CFLAGS=-std=c++0x -O3 -march=native -mtune=native -Wall -c -fmessage-length=0 -mpopcnt
LDFLAGS= -lpng12 -lm -lz -lpthread
SOURCES=Source.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=DM

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) -pg -o $@ $(OBJECTS) $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@



#target: dependencies
#[tab] system command
