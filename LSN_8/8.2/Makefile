CC = g++
CFLAGS = -Wall -O3 -std=c++11
LIBS = -larmadillo

main.exe : main.o random.o lib.o functions.o annealing.o
	$(CC) random.o lib.o functions.o annealing.o main.o -o main.exe $(LIBS) 

main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

lib.o : lib.cpp lib.h
	$(CC) -c lib.cpp -o lib.o $(CFLAGS)

functions.o : functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o $(CFLAGS)

annealing.o : annealing.cpp annealing.h
	$(CC) -c annealing.cpp -o annealing.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out

run : main.exe
	./main.exe
