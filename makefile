all: main_iplib.c ip_lib/ip_lib.o bmp/bmp.o
	gcc ip_lib/ip_lib.o bmp/bmp.o main_iplib.c -o image_editor -Wall -lm

ip_lib/ip_lib.o: ip_lib/ip_lib.c bmp/bmp.o
	gcc ip_lib/ip_lib.c bmp/bmp.o -o ip_lib/ip_lib.o -Wall -lm -c

bmp/bmp.o: bmp/bmp.c
	gcc bmp/bmp.c -o bmp/bmp.o -Wall -lm -c

clean: 
	rm image_editor ip_lib/ip_lib.o bmp/bmp.o
