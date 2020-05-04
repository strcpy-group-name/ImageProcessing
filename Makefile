CC=gcc
CFLAGS= --ansi --pedantic -Wall -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra
LIBFLAGS= -Wall

main: main.o bmp.o bmp.h
	$(CC) main.o bmp.o bmp.h $(CFLAGS) -o$@

main.o: test_bmp.c bmp.h
	$(CC) test_bmp.c $(LIBFLAGS) -c -o$@

bmp.o: bmp.c bmp.h
	$(CC) bmp.c $(LIBFLAGS) -c -o$@

clean:
		@rm -f main
		@rm -f main.o
		@rm -f bmp.o
