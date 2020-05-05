CC=gcc
CFLAGS= --ansi --pedantic -Wall -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra -ggdb -g
LIBFLAGS= -Wall -ggdb -g

test_ipmat: test_ipmat.o ip_lib.o ip_lib.h bmp.o bmp.h
	$(CC) test_ipmat.o ip_lib.o bmp.o $(CFLAGS) -o$@

test_ipmat.o: test_ipmat.c ip_lib.h bmp.h
	$(CC) test_ipmat.c $(LIBFLAGS) -c -o$@

bmp.o: bmp.c bmp.h
	$(CC) bmp.c $(LIBFLAGS) -c -o$@

ip_lib.o: ip_lib.c ip_lib.h
	$(CC) ip_lib.c $(LIBFLAGS) -c -o$@

clean:
		@rm -f bmp.o
		@rm -f test_ipmat
		@rm -f test_ipmat.o
		@rm -f ip_lib.o
