CC = g++
CFLAGS = -lm -lgsl -lblas
OBJ = Main.o System.o
DEPS = $(OBJ) System.h
DBG = -g

Main: $(DEPS)
	$(CC) -o $@ $(OBJ) $(CFLAGS)

debug: $(DEPS)
	$(CC) -o Main $(OBJ) $(CFLAGS) $(DBG)

plot: Plot.gp
	gnuplot Plot.gp -p

clean:
	rm *.o Main
