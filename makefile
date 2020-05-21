all: src/main_iplib.c build/ip_lib/ip_lib.o build/bmp/bmp.o
	gcc build/ip_lib/ip_lib.o build/bmp/bmp.o src/main_iplib.c -o build/image_editor -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra

build/ip_lib/ip_lib.o: src/ip_lib/ip_lib.c build/bmp/bmp.o build/ip_lib
	gcc src/ip_lib/ip_lib.c -o build/ip_lib/ip_lib.o -c -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra

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
