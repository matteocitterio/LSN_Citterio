CC = mpicxx
CFLAGS = -Wall -O3 -std=c++11
LIBS = -larmadillo

main.exe : main.o random.o geneticalgorithm.o
	$(CC) random.o geneticalgorithm.o main.o -o main.exe $(LIBS) 

main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

geneticalgorithm.o : geneticalgorithm.cpp geneticalgorithm.h
	$(CC) -c geneticalgorithm.cpp -o geneticalgorithm.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out

run : main.exe
	./main.exe
