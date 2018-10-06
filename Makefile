CC = g++
CFLAGS = -Wall

pso: pso.cpp
	 $(CC) $(CFLAGS) -o $@ pso.cpp

clean:
	 rm -f pso