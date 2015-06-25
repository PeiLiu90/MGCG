CXX=g++

objects = main.o

run: $(objects)
	g++ -o run $(objects)

.PHONY:clean
clean:
	-rm run $(objects)
