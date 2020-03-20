
all:
	gcc -Wall -O3 -g -mtune=native -march=native -fopenmp oplib.c 
clean:
	rm -f a.out
