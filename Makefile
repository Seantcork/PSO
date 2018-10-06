CC = g++
CFLAGS = -Wall

pso: pso.cpp
	 $(CC) $(CFLAGS) -o $@ pso

clean:
	 rm -f pso