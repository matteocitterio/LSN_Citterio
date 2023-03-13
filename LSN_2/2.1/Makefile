CC = g++
CFLAGS = -Wall -O3 --std=c++11 -I/opt/homebrew/Cellar/armadillo/12.0.1/include
CFLAGSNOARM = -Wall -O3 --std=c++11
LDFLAGS = -L/opt/homebrew/Cellar/armadillo/12.0.1/lib
LIBS = -larmadillo

main.exe : main.o random.o functions.o lib.o integral.o
	$(CC) random.o functions.o lib.o integral.o main.o -o main.exe $(CFLAGSNOARM)
main.o : main.cpp
	$(CC) -c main.cpp -o main.o 
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o 
functions.o: functions.cpp functions.h
	$(CC) -c functions.cpp -o functions.o 
lib.o : lib.cpp lib.h
	$(CC) -c lib.cpp -o lib.o 
integral.o : integral.cpp integral.h
	$(CC) -c integral.cpp -o integral.o 
clean :
	rm *.o main.exe seed.out
run : main.exe
	./main.exe




