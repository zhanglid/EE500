CC = gcc
LN = gcc
All : prob1

ASAM2D:
	@g++ -c ASAM-Main.cpp
	@g++ -o ASAM-Main ASAM-Main.o
	./ASAM-Main

clean:
	rm -f *.o ASAM-Main
	rm -rf Cauchy Gauss Sinc
	rm -f InputFxn.dat Errors.dat
