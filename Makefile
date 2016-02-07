CC = g++
CFLAGS = -lm -lgsl -lblas

Main: Main.o System.o
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm *.o Main
