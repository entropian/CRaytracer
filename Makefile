CC=g++
LFLAGS =-lGL -lGLEW -lglfw -lGLU
CFLAGS = -std=c++11 -c

all:
		$(CC) -o craytracer main.cpp $(LFLAGS)
