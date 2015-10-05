# DE_Self-Consistent
Differential Evolution for self-consistent equations

# DE_Self-Consistent (Version 0.9.0)
=========================
Author: Ehsan Zahedinejad ezahedin@ucalgary.ca
=========================

Title:
=========================
Differential Evolution for solving self-consistent equation.
Program can be ran on single and mutiple processors.
User should define the mesh size as well as the function that DE should evaluate at each iteration.

Libraries:
============================================
MPI (Message Passing Interface), MKL (Math Kernal Library)
A single thread compilation of mkl leads to better efficiency. This is why I compiled using a single thread command.

Compilation
======================
Compilation command: ./make
However you should alway make sure that the program refers to the right location of libraries.
Edit makefile according to the location of you libraries.


