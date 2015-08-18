all:
	gcc -g -o main main.c -lm -Wall

run:
	./main

clean:
	rm main
