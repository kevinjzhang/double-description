all: main.cc
	g++ -std=c++17 -O3 `regina-engine-config --cflags --libs` main.cc -o double-description
	
