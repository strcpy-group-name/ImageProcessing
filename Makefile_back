CC=gcc
CFLAGS= -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra
LIBFLAGS= -Wall -lm
DBGFLAGS= -O1 -g -Wall -lm

all: main test debug

bmp.o: bmp.c bmp.h
	$(CC) bmp.c $(LIBFLAGS) -c -o$@

ip_lib.o: ip_lib.c ip_lib.h
	$(CC) ip_lib.c $(CFLAGS) -c -o$@

main: main_iplib.o ip_lib.o ip_lib.h bmp.o bmp.h
	$(CC) main_iplib.o ip_lib.o bmp.o $(CFLAGS) -o$@

main_iplib.o: main_iplib.c ip_lib.h bmp.h
	$(CC) main_iplib.c $(CFLAGS) -c -o$@

test: test.o ip_lib.o ip_lib.h bmp.h bmp.o
	$(CC) test.o ip_lib.o bmp.o $(CFLAGS) -o$@

test.o: test.c ip_lib.h bmp.h
	$(CC) test.c $(CFLAGS) -c -o$@

debug.o: test.c ip_lib.h bmp.h
	$(CC) test.c $(DBGFLAGS) -c -o$@

debug: debug.o ip_lib.h ip_lib.o bmp.h bmp.o
	$(CC) debug.o ip_lib.o bmp.o $(CFLAGS) -o$@

clean:
		@rm -f test.o
		@rm -f test
		@rm -f bmp.o
		@rm -f main
		@rm -f main.o
		@rm -f ip_lib.o
		@rm -f main_iplib.o
		@rm -f debug
		@rm -f debug.o
