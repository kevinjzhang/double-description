all: main.cc
	g++ -O3 `regina-engine-config --cflags --libs` main.cc -o double-description
	
