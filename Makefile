CC=gcc
CFLAGS= -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra
LIBFLAGS= -Wall -lm

main_iplib: main_iplib.o ip_lib.o ip_lib.h bmp.o bmp.h
	$(CC) main_iplib.o ip_lib.o bmp.o $(CFLAGS) -o main_iplib

bmp.o: bmp.c bmp.h
	$(CC) bmp.c $(LIBFLAGS) -c -o bmp.o

ip_lib.o: ip_lib.c ip_lib.h
	$(CC) ip_lib.c $(CFLAGS) -c -o ip_lib.o

main_iplib.o: main_iplib.c ip_lib.h bmp.h
	$(CC) main_iplib.c $(CFLAGS) -c -o main_iplib.o

clean:
		@rm -f bmp.o
		@rm -f main_iplib
		@rm -f ip_lib.o
		@rm -f main_iplib.o