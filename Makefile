all: main.cc
	g++ -fopenmp -O3 -std=c++17 `regina-engine-config --cflags --libs` main.cc -o double-description
	
