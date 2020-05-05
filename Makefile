CC=gcc
CFLAGS= --ansi --pedantic -Wall -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra
LIBFLAGS= -Wall

main_test_ipmat: test_ipmat.o ip_lib.o ip_lib.h
	$(CC) test_ipmat.o ip_lib.o ip_lib.h $(CFLAGS) -o$@

main: main.o bmp.o bmp.h
	$(CC) main.o bmp.o bmp.h $(CFLAGS) -o$@

test_ipmat.o: test_ipmat.c ip_lib.h
	$(CC) test_ipmat.c $(LIBFLAGS) -c -o$@

main.o: test_bmp.c bmp.h
	$(CC) test_bmp.c $(LIBFLAGS) -c -o$@

bmp.o: bmp.c bmp.h
	$(CC) bmp.c $(LIBFLAGS) -c -o$@





ip_lib.o: ip_lib.c ip_lib.h
	$(CC) ip_lib.c $(LIBFLAGS) -o$@

clean:
		@rm -f main
		@rm -f main.o
		@rm -f bmp.o
		@rm -f test_ipmat
		@rm -f test_ipmat.O
		@rm -f ip_lib.o
