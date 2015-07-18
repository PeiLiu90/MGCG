objects = main.o
CXX=g++
CFLAGS=-O3 

run: $(objects)
	$(CXX) $(CFLAGS) -o run $(objects)

main.o: main.cpp MGCG.h MultiGrid.h RBSGS.h vector2d.h vec.h
	$(CXX) $(CFLAGS) -c main.cpp -o main.o
.PHONY:clean
clean:
	-rm run $(objects)
