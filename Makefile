LIBS = /usr/local/lib64/libregina-engine.so.5.1
LD_LIBRARY_PATH=/usr/local/lib64/

export LD_LIBRARY_PATH=/usr/local/lib64/

all: main.cc
	g++ `regina-engine-config --cflags --libs` main.cc -o double-description
	
