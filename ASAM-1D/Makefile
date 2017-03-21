CC = gcc
LN = gcc
All : prob1

ASAM1D:
	@g++ -c SAMUtils.cpp
	@g++ -c Main.cpp
	@g++ -o Main Main.o SAMUtils.o
	./Main

clean:
	rm -f *.o Main SAMUtils
	rm -rf Cauchy Laplace Gauss Sinc Tanh Tri
	rm -f InputFxn.dat Errors.dat
