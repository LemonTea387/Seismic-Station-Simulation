all: build clean main

main:
	mpicc src/main.c src/seismic_sensor.c src/util.c src/struct_data.c src/base.c -Wall -fopenmp -lm -o bin/app
	mpirun -np 10 --oversubscribe ./bin/app 3 3 5

clean:
	rm bin/*

build:
	mkdir -p ./bin