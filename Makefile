CC=gcc
CFLAGS= --ansi --pedantic -Wall -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra -ggdb -g
LIBFLAGS= -Wall -ggdb -g

main.out: main.o ip_lib.o ip_lib.h bmp.o bmp.h
	$(CC) main.o ip_lib.o bmp.o $(CFLAGS) -o$@

main.o: main.c ip_lib.h bmp.h
	$(CC) main.c $(LIBFLAGS) -c -o$@

bmp.o: bmp.c bmp.h
	$(CC) bmp.c $(LIBFLAGS) -c -o$@

ip_lib.o: ip_lib.c ip_lib.h
	$(CC) ip_lib.c $(LIBFLAGS) -c -o$@

clean:
		@rm -f bmp.o
		@rm -f main
		@rm -f main.o
		@rm -f ip_lib.o
