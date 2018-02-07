CC = gcc -g
CFLAGS = `sisl-config --cflags` 
LIBS = ../src/.libs/libsisl.a `sisl-config --libs` `wmpi-config --libs` `gsl-config --libs`

test: test.c
	$(CC) $(CFLAGS) -o test test.c $(LIBS)
