
all:
	cd src;	gcc -Wall -O3 -g  -mtune=native -march=native -ffast-math -fopenmp -funroll-loops oplib.c

amd:
	cd src; clang -Wall -g -O3 -mtune=znver2 -march=znver2 -ffast-math -funroll-loops  -fopenmp oplib.c

clean:
	rm -f a.out
