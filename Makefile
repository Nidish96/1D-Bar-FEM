CC = g++
CFLAGS = -lm -lgsl -lblas
OBJ = Main.o System.o
DEPS = $(OBJ) System.h
DBG = -g

Nodes = 20
Elements = -1
Interval = 500

Main: $(DEPS)
	$(CC) -o $@ $(OBJ) $(CFLAGS)

debug: $(DEPS)
	$(CC) -o Main $(OBJ) $(CFLAGS) $(DBG)

plot: Plot.gp
	./Main $(Nodes) $(Elements) $(Interval)>OUTPUT.dat && gnuplot5 Plot.gp -p

clean:
	rm *.o Main
