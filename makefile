CC = mpic++# Compiler to use
OPTIONS = -O3 -Wall -std=c++0x -pthread -fopenmp -fpermissive# -g for debug, -O2 for optimise and -Wall additional messages
INCLUDES = -I . # Directory for header file
OBJS = GNG_Main.o GrowingNetwork.o LinedList.o Neuron.o ParameterParser.o # List of objects to be build
#.PHONY: all clean # To declare all, clean are not files
 
all: $(OBJS)
	@echo "Building.." # To print "Building.." message
	$(CC) $(OPTIONS) $(INCLUDES) $(OBJS) -o ./run/gng_linux 
 
%.o: ./src/%.cpp    
	$(CC) $(OPTIONS) -c ./src/$*.cpp $(INCLUDES) 
#	$(CC) $(OPTIONS) -c $< -o $@
list:
	@echo $(shell ls) # To print output of command 'ls'
 
clean:
	@echo "Cleaning up.."
	-rm -rf *.o # - prefix for ignoring errors and continue execution
	-rm ./run/gng_linux

