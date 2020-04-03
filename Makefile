
all:
	#gcc -Wall -O3 -g  -mtune=native -march=native -fopenmp oplib.c
	clang -O3 -mtune=znver2 -march=znver2 -ffast-math -funroll-loops  -fopenmp test_prll.c
clean:
	rm -f a.out
