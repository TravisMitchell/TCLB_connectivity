
all : main

main: main.o lodepng.o
	g++ -o $@ $^

%.o : %.cpp
	g++ -O3 -c -o $@ $<