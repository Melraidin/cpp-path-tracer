path-tracer: main.cpp ray.h sphere.h sphere.cpp enums.h
	g++ -O3 -fopenmp main.cpp sphere.cpp -I /usr/include/eigen3 -lpng -o path-tracer

# sphere: sphere.h sphere.cpp
# 	g++ sphere.cpp -o

eigen-test:
	g++ eigen-test.cpp -I /usr/include/eigen3 -o eigen-test

png-write-test:
	g++ png-write-test.cpp -lpng -o png-write-test

smallpt: smallpt.cpp
	g++ -O3 -fopenmp smallpt.cpp -o smallpt
