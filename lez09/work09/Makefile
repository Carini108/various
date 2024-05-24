CC = g++
CFLAGS = -Wall -O3 --std=c++11

GA.exe : GA.o random.o
	$(CC) random.o GA.o -o GA.exe
GA.o : GA.cpp
	$(CC) -c GA.cpp -o GA.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o GA.exe
cleanoutput :
	rm *.genetic
