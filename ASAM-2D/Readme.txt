Author (version 2.1) : Olaoluwa Adigun

Trains 6 SAMs to approximate a user-specified 2-dimensional function using samples of the objective function.
The additional features are:
--- Computing the conditonal variance for the Fuzzy approximation.
--- Fuzzy aproximation for noisy functions

How to run:
===========
1) - The program takes a function f(x) and finds its fuzzy approximation. You can specify the desired function.
x	y	f(x,y)
for example:
5.56565657	6.14141414	0.0929060987	

2) - The program outputs 3 folders: Gauss, Sinc, and Cauchy. These contain the ASAM approximation details for the three different set functions. You should tune the number of rules, number of iterations, learning rates, and initializations to get better performance on your approximations. It also outputs the "Errors.dat" file. This is a log of MSEs for all the ASAMs at different points in the adaptation. "InputFxn.dat" is just a sanity-check procedure. 

3) There is a make file for running the code. 
----- From the terminal or command prompt, navigate to the location of the software (ASAM-2D).  
-----Type "make ASAM2D" on the CMD terminal. 
-----Use "make clean" to delete all the output files (Optional).

Compilation Requirements: (Windows) 
===================================
1) Eigen library: You need to download the Eigen matix manipulation library from http://eigen.tuxfamily.org/dox/
Your compiler needs to know where you put the library. You do this in Visual Studio by changing the "project properties->VC++ Directories->Include Directories". *Add* the eigen library folder path to the include directories.

Caveats:
========
1) - This 2D ASAM library differs from the previous version. I am using large arrays instead of matrices to store and manipulate f(x,y). You may need to do some further heavy lifting if you want to train on >1000 samples. 
2) - The program erases outputs from previous runs every time you run it.
