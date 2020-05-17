#all: src/main_iplib.c build/ip_lib/ip_lib.o build/bmp/bmp.o
#	gcc build/ip_lib/ip_lib.o build/bmp/bmp.o src/main_iplib.c -o build/image_editor -Wall -lm

all: src/test_main.c build/ip_lib/ip_lib.o build/bmp/bmp.o
	gcc build/ip_lib/ip_lib.o build/bmp/bmp.o src/test_main.c -o build/testrun -Wall -lm

debug: src/test_main.c build/ip_lib/ip_lib.o buildbmp/bmp.o
	gcc build/ip_lib/ip_lib.o build/bmp/bmp.o src/test_main.c -o build/debugrun -ansi -pedantic -Wall -lm -g

build/ip_lib/ip_lib.o: src/ip_lib/ip_lib.c build/bmp/bmp.o build/ip_lib
	gcc src/ip_lib/ip_lib.c build/bmp/bmp.o -o build/ip_lib/ip_lib.o -ansi -pedantic -Wall -lm -c

build/ip_lib: build
	mkdir -p build/ip_lib

build/bmp/bmp.o: src/bmp/bmp.c build/bmp
	gcc src/bmp/bmp.c -o build/bmp/bmp.o -Wall -lm -c

build/bmp: build
	mkdir -p build/bmp

build: 
	mkdir -p build

clean: 
	rm -Rf build
