
all:
	gcc -Wall -O2 -g  -mtune=native -march=native -ffast-math -fopenmp -funroll-loops oplib.c

amd:
	clang -Wall -g -O2 -mtune=znver2 -march=znver2 -ffast-math -funroll-loops  -fopenmp oplib.c

clean:
	rm -f a.out
