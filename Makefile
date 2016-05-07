CC=g++
LFLAGS =-lGL -lGLEW -lglfw -lGLU
CFLAGS = -std=c++11 -c

all:
		$(CC) -g -o craytracer main.cpp $(LFLAGS)
