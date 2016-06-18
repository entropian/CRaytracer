CC=g++
LFLAGS = -lGL -lGLEW -lglfw -lGLU
CFLAGS = -std=c++11 -c

OS = $(shell uname)

ifeq ($(OS), Darwin)
all:
		$(CC) -g -o craytracer main.cpp -lGLEW -lglfw3 -framework Cocoa -framework OpenGL -framework IOKit -framework CoreVideo
else
all:
		$(CC) -g -o craytracer main.cpp $(LFLAGS)
endif
