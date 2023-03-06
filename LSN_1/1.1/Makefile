CC = g++
CFLAGS = -Wall -O3 --std=c++11 -I/opt/homebrew/Cellar/armadillo/12.0.1/include
LDFLAGS = -L/opt/homebrew/Cellar/armadillo/12.0.1/lib
LIBS = -larmadillo

main.exe : main.o random.o funzioni.o
	$(CC) random.o funzioni.o main.o -o main.exe $(LIBS) $(LDFLAGS)
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
run : main.exe
	./main.exe




