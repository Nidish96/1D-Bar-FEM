CC = g++
CFLAGS = -lm -lgsl -lblas

Main: Main.o System.o
	$(CC) -o $@ $^ $(CFLAGS)

debug: Main.o System.o
	$(CC) -o Main $^ $(CFLAGS) -g

clean:
	rm *.o Main
