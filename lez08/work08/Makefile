CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : ./prng/random.cpp ./prng/random.h stats.h
	$(CC) -c ./prng/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
cleandat :
	rm *.dat