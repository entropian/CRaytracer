CC=g++
LFLAGS = -lGL -lGLEW -lglfw -lGLU -lpthread
CFLAGS = -std=c++11 -c

OS = $(shell uname)

ifeq ($(OS), Darwin)
all:
		$(CC) -g -o craytracer main.cpp -lGLEW -lglfw3 -framework Cocoa -framework OpenGL -framework IOKit -framework CoreVideo
else
all:
		$(CC) -g -o craytracer main.cpp aabb.cpp camera.cpp lights.cpp materials.cpp mesh.cpp noise.cpp sampling.cpp texture.cpp shapes/box.cpp shapes/cylinder.cpp shapes/disk.cpp shapes/generic.cpp shapes/instanced.cpp shapes/objecttype.cpp shapes/plane.cpp shapes/rect.cpp shapes/shapes.cpp shapes/sphere.cpp shapes/torus.cpp shapes/triangle.cpp util/ray.cpp util/math.cpp util/util.cpp $(LFLAGS)
endif
