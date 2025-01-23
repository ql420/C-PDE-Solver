IN = 1
T1 = --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 1 --in $(IN)
T2 = --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 2 --in $(IN)
T3 = --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 3 --in $(IN)
T4 = --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 4 --in $(IN)
TARGET = myprog
CC = g++
CXXFLAGS = -Wall -g -O2 -c
LDLIBS = -lblas -fopenmp -lboost_program_options
default: 

	$(CC) $(CXXFLAGS) -o main.o main.cpp
	$(CC) $(CXXFLAGS) -o ShallowWater.o ShallowWater.cpp
	$(CC) main.o ShallowWater.o $(LDLIBS) -o myprog

test1: $(TARGET)
	./$(TARGET) $(T1)

test2: $(TARGET)
	./$(TARGET) $(T2)

test3: $(TARGET)
	./$(TARGET) $(T3)

test4: $(TARGET)
	./$(TARGET) $(T4)

.PHONY:clean

clean:
	-rm -f *.o
