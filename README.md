# fastQSC_twistedHagedorn
A program for computing the Hagedorn temperature in N=4 SYM in the presence of a chemical potential for one of the R-symmetries.

This is a sample of the code used in the preparation of "Twisting the Hagedorn temperature in planar N=4 super Yang-Mills" by J.Minahan, S. Ekhammar and C.Thull, available at https://arxiv.org/abs/2512.05810.

This code is written in C++. Please use your favourite compiler to produce an executable for your system. For the compilation of this code the Class Library for Numbers is required, see https://www.ginac.de/CLN/ for download and documentation. On Linux the compilation can be done with g++ using the commmand "g++ auxmathfuncs.cpp interpolation.cpp linSolve.cpp QSC.cpp main.cpp -lcln -lm -o test.out".

Once compiled, the program can be executed from the command line. On Linux the corresponding command would be "./test.out". The execution will produce 3 text files:
1) messages_omega0.707107.txt informing on the progress of the execution, giving the error norm and the value of T_H at each step of the optimization
2) rangeCoupling_omega0.707107.txt containing the coupling values and corresponding data arrays at full precision as found by the algorithm
3) rangeCoupling_omega0.707107_onlyT_double.txt containing only the coupling values and exp(-1/(4*T_H)) with double precision

If a number n of values has already been computed at a previous execution, these can be loaded at a repeat execution: copy the two datafiles from the previous execution into the "Old" subfolder and then run the program with the command "./test.out n " where n is the number of values already computed.
