objects = main.o
CXX=g++
CFLAGS=-O3 

run: $(objects)
	$(CXX) $(CFLAGS) -o run $(objects)

main.o: main.cpp MGCG3d.h MultiGrid3d.h RBSGS3d.h vector3d.h vec.h
	$(CXX) $(CFLAGS) -c main.cpp -o main.o
.PHONY:clean
clean:
	-rm run $(objects)
