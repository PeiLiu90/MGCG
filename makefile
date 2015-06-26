CXX=g++
CFLAGS=-O3 

objects = main.o

all: $(objects) run

$(objects): %.o: %.cpp
	$(CXX) -c $(CFLAGS) $< -o $@ 

run: $(objects)
	$(CXX) $(CFLAGS) -o run $(objects)

.PHONY:clean
clean:
	-rm -f run $(objects)
