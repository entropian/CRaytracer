CC=g++
LFLAGS =-lGL -lGLEW -lglfw -lGLU
CFLAGS = -std=c++11 -c

all:
		$(CC) -o craytracer  util/vec.h util/mat.h aabb.h buildscene.h camera.h lights.h materials.h  sampling.h  sceneobj.h  shading.h gl/glcode.h gl/glshaders.h shapes/aabox.h shapes/cylinder.h shapes/disk.h shapes/plane.h shapes/rect.h shapes/shapes.h shapes/sphere.h shapes/torus.h shapes/triangle.h util/constants.h util/math.h util/ray.h util/shaderec.h util/vec.h main.cpp $(LFLAGS)
