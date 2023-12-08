# LDLtr
                                                                       
 LDLtr: An LDLT trust-region region algorithm for the minimization 
 of a general nonlinear function when only gradients are available.                                     
                                                                       
                            minimize f(x)                                   
 
 The methods seek a point "norm(grad(x)) = 0", where grad(x) is
 the gradient of the objective f(x).
 
## Setup
Files are organized in the following way:

	FORTRAN/
 		cutest/			(files for integration with CUTEst)
 		examples/		(stand-alone example)
   		funcs/			(objective function(s) for example)     		
   		include/		(source: triang. solves, LDL updates, and more)
     	main/			(source: main algorithm, parameters)
       	objs/ 			(complied files with ending: *.o and *.mod)
	 	test/ 			(internal tests)

	MATLAB/
 		Statistics/		(output of algorithm runs)
   		EXAMPLE.m		(stand-alone example)
     	LDLtr.m			(main algorithm)
       	runUCset.m		(experiment with CUTEst problems)
		more files
       		
In order to run the experiments based on the large set of problems the
testing library CUTEst (https://github.com/ralna/CUTEst/wiki) has to be installed. 

### Fortran
The ldl algorithm can be re-complied from the source files or one can make use
of the precomplied binanries. For the simplest experience on a unix machine,
or on Windows (with Cygwin, or similar) one can first navigate to
FORTRAN/examples and then run

		./test_ldlqn

In order to interface the solver with CUTEst (after CUTEst has been installed)
one can copy the files within the FORTRAN/cutest folder to the corresponding
folders of the CUTEst installation. For instance, one can copy the files
from FORTRAN/cutest/src/ldltr to "LOCAL_CUTEST"/cutest/scr/ here 
"LOCAL_CUEST" refers to the local installation folder. When LDLtr is interfaced 
with CUTEst one can run the algorithm with e.g.,

		runcutest -p ldltr -D BDQRTIC

where BDQRTIC is a problem name from CUTEst. This assumes (as is usually the case)
that the master list with all CUTEst problems is on the search path. 

To compile a program with the algorithm one will have to link to some of the
solver libraries. For instance, to recompile the program with the Fortran
example (e.g., after making changes to the size of problem "n" in 
FORTRAN/examples/test_ldlqn.f90), one first compiles the program

		 gfortran -c test_ldlqn.f90 -I ../objs/ -J ../objs/

Second, one can link it to all the relevant algorithm libraries

		gfortran -o test_ldlqn test_ldlqn.o ../objs/kind_parameter.o ../objs/ldlqn_bfgsRDRup.o ../objs/ldlqn_invLDLCG.o ../objs/ldlqn.o ../objs/ldlqn_pars.o ../objs/ldlqn_triang.o ../objs/linesearch.o ../objs/rosenbrock.o

For reference, the relevant algorithm libraries are in FORTRAN/objs/
and are: kind_parameter.o, ldlqn_bfgsRDRup.o, ldlqn_invLDLCG.o, 
ldlqn.o, ldlqn_pars.o, ldlqn_triang.o, and linesearch.o

An example Rosenbrock objective function is in objs/rosenbrock.o 

After compiling and linking the program you can run it by ```./test_ldlqn ```
		
### Matlab
With CUTEst for Matlab installed one can modify the path specified
in MATLAB/runUCset.m to point ot the local CUTEst installation.

## Examples

### Fortran
To try-out an example you can navigate to FORTRAN/examples and call

		./test_ldlqn

The result on linux computer with gfortran 10.5 is

```
************************************************************************
 *
 * LDLTR1
 * 
 * LDL Trust-Region Quasi-Newton Algorithm
 *
 * n=    50
 * fx=25.000
 * gx=10.000
 *
 ************************************************************************
 
itn     obj       ng        delt      ns        rho    trit   phi    dn   
     1 0.250E+02 0.100E+02 0.200E+01 0.100E+01   0.00      0   0.00   8.53
     2 0.112E+02 0.296E+02 0.400E+01 0.200E+01   0.90      3  -0.00   2.53
     3 0.112E+02 0.296E+02 0.200E+01 0.562E+00  -0.02      4  -0.00   2.53
     4 0.899E+01 0.172E+02 0.200E+01 0.268E+00   0.64      4  -0.00   1.73
     5 0.899E+01 0.172E+02 0.100E+01 0.186E+01 -50.97      3  -0.00   1.73
     6 0.899E+01 0.172E+02 0.500E+00 0.542E+00  -0.59      5  -0.00   1.73
     7 0.894E+01 0.435E+02 0.250E+00 0.500E+00   0.03      5  -0.00   1.72
     8 0.894E+01 0.435E+02 0.125E+00 0.181E+00  -0.01      4  -0.00   1.72
     9 0.698E+01 0.191E+02 0.125E+00 0.112E+00   0.63      3  -0.00   1.19
    10 0.628E+01 0.103E+02 0.125E+00 0.125E+00   0.61      3  -0.00   1.70
    11 0.600E+01 0.145E+02 0.125E+00 0.125E+00   0.69      3  -0.00   2.09
    12 0.568E+01 0.163E+02 0.250E+00 0.125E+00   0.80      3  -0.00   2.00
    13 0.539E+01 0.250E+02 0.250E+00 0.250E+00   0.36      3   0.00   2.06
    14 0.500E+01 0.271E+02 0.250E+00 0.250E+00   0.40      4   0.00   2.51
itn     obj       ng        delt      ns        rho    trit   phi    dn   
    15 0.500E+01 0.271E+02 0.125E+00 0.250E+00  -0.51      3  -0.00   2.51
    16 0.466E+01 0.207E+02 0.125E+00 0.845E-01   0.25      3  -0.00   2.34
    17 0.434E+01 0.199E+02 0.125E+00 0.125E+00   0.39      4  -0.00   2.49
    18 0.392E+01 0.171E+02 0.125E+00 0.125E+00   0.54      3  -0.00   2.57
    19 0.356E+01 0.131E+02 0.125E+00 0.125E+00   0.68      3  -0.00   2.60
    20 0.335E+01 0.142E+02 0.125E+00 0.125E+00   0.59      3   0.00   2.59
    	.		.		.		.		.
	.		.		.		.		.
	.		.		.		.		.
	.		.		.		.		.	
itn     obj       ng        delt      ns        rho    trit   phi    dn   
   120 0.717E-06 0.170E-01 0.625E-01 0.174E-02   1.24      0  -0.00   0.77
   121 0.474E-06 0.536E-02 0.625E-01 0.657E-03   1.09      0  -0.00   0.79
   122 0.438E-06 0.419E-02 0.625E-01 0.231E-03   0.93      0  -0.00   0.77
   123 0.412E-06 0.306E-02 0.625E-01 0.100E-03   1.57      0  -0.00   0.70
   124 0.337E-06 0.450E-02 0.625E-01 0.256E-03   1.53      0  -0.00   0.68
   125 0.206E-06 0.800E-02 0.625E-01 0.458E-03   1.41      0  -0.00   0.69
   126 0.899E-07 0.689E-02 0.625E-01 0.536E-03   1.28      0  -0.00   0.73
   127 0.444E-07 0.277E-02 0.625E-01 0.324E-03   1.13      0  -0.00   0.71
   128 0.368E-07 0.861E-03 0.625E-01 0.991E-04   1.10      0  -0.00   0.67
   129 0.347E-07 0.671E-03 0.625E-01 0.370E-04   1.58      0  -0.00   0.62
   130 0.295E-07 0.125E-02 0.625E-01 0.728E-04   1.57      0  -0.00   0.60
   131 0.194E-07 0.201E-02 0.625E-01 0.120E-03   1.50      0  -0.00   0.62
   132 0.815E-08 0.192E-02 0.625E-01 0.156E-03   1.39      0  -0.00   0.66
   133 0.228E-08 0.918E-03 0.625E-01 0.121E-03   1.26      0  -0.00   0.66
   134 0.117E-08 0.178E-03 0.625E-01 0.492E-04   1.17      0  -0.00   0.62
itn     obj       ng        delt      ns        rho    trit   phi    dn   
   135 0.104E-08 0.112E-03 0.625E-01 0.151E-04   1.33      0  -0.00   0.61
   136 0.958E-09 0.151E-03 0.625E-01 0.863E-05   1.65      0  -0.00   0.64
   137 0.725E-09 0.240E-03 0.625E-01 0.184E-04   1.52      0  -0.00   0.70
   138 0.410E-09 0.278E-03 0.625E-01 0.236E-04   1.47      0  -0.00   0.72
   139 0.140E-09 0.186E-03 0.625E-01 0.257E-04   1.33      0  -0.00   0.66
   140 0.516E-10 0.595E-04 0.625E-01 0.155E-04   1.22      0  -0.00   0.61
 
Status:  0

```


With CUTEst installed one can issue a command such as 

		jburst@francis:~/Dropbox/git/LDLtr/FORTRAN$ runcutest -p ldltr -D BDQRTIC

The output on linux computer with gfortran 10.5 is

```
jburst@francis:~/Dropbox/git/LDLtr/FORTRAN$ runcutest -p ldltr -D BDQRTIC
sifdecoder -A pc64.lnx.gfo -st   BDQRTIC

 Problem name: BDQRTIC

 Double precision version will be formed

 The objective function uses 4996 linear groups
 The objective function uses 4996 nonlinear groups
 
 There are 5000 free variables
 
 
 File successfully decoded
 CUTEST: tools (double precision version) compiled successfully
 CUTEst: ldltr (double precision version) compiled successfully
 ************************************************************************
 *
 * LDLTR1
 * 
 * LDL Trust-Region Quasi-Newton Algorithm
 *
 * n=  5000
 * fx=******
 * gx=******
 *
 ************************************************************************
 
itn     obj       ng        delt      ns        rho    trit   phi    dn   
     1 0.113E+07 0.150E+07 0.213E+01 0.106E+01   0.00      0   0.00 100.00
     2 0.100E+06 0.905E+05 0.213E+01 0.246E+02   0.71      2  -0.00 100.00
     3 0.508E+05 0.141E+06 0.213E+01 0.171E+02   0.60      3  -0.00 100.00
     4 0.336E+05 0.782E+05 0.425E+01 0.189E+01   1.53      0  -0.00 100.00
     5 0.232E+05 0.303E+05 0.425E+01 0.282E+01   1.31      0  -0.00 100.00
     6 0.208E+05 0.133E+05 0.425E+01 0.120E+01   1.37      0  -0.00 100.01
     7 0.201E+05 0.420E+04 0.425E+01 0.453E+00   1.27      0  -0.00 100.00
     8 0.200E+05 0.501E+03 0.425E+01 0.182E+00   1.10      0  -0.00  99.99
     9 0.200E+05 0.327E+02 0.425E+01 0.119E+00   0.85      0  -0.00  99.74
    10 0.200E+05 0.400E+02 0.425E+01 0.504E-01   1.09      0  -0.00  97.54
    11 0.200E+05 0.589E+01 0.425E+01 0.166E-01   1.31      0  -0.00  82.95
    12 0.200E+05 0.370E+01 0.425E+01 0.155E-01   1.17      0  -0.00  57.20
    13 0.200E+05 0.306E+01 0.425E+01 0.100E-01   1.64      0  -0.00  26.02
    14 0.200E+05 0.159E+01 0.425E+01 0.285E-01   1.29      0  -0.00  23.37
itn     obj       ng        delt      ns        rho    trit   phi    dn   
    15 0.200E+05 0.344E+00 0.425E+01 0.144E-01   1.29      0  -0.00  24.60
    16 0.200E+05 0.349E+00 0.425E+01 0.622E-02   1.36      0  -0.00  26.14
    17 0.200E+05 0.442E+00 0.425E+01 0.527E-02   1.49      0  -0.00  19.61
    18 0.200E+05 0.383E+00 0.425E+01 0.751E-02   1.37      0  -0.00  14.13
    19 0.200E+05 0.170E+00 0.425E+01 0.549E-02   1.34      0  -0.00  11.10
    20 0.200E+05 0.969E-01 0.425E+01 0.323E-02   1.42      0  -0.00   8.11
    21 0.200E+05 0.167E+00 0.425E+01 0.342E-02   1.54      0  -0.00   4.93
    22 0.200E+05 0.232E+00 0.425E+01 0.633E-02   1.50      0  -0.00   3.12
    23 0.200E+05 0.206E+00 0.425E+01 0.879E-02   1.42      0  -0.00   2.47
    24 0.200E+05 0.964E-01 0.425E+01 0.719E-02   1.31      0  -0.00   2.55
    25 0.200E+05 0.358E-01 0.425E+01 0.233E-02   1.23      0  -0.00   2.74
    26 0.200E+05 0.249E-01 0.425E+01 0.642E-03   1.19      0  -0.00   2.69
    27 0.200E+05 0.110E-01 0.425E+01 0.326E-03   1.13      0  -0.00   2.74
    28 0.200E+05 0.200E-02 0.425E+01 0.529E-04   1.06      0  -0.00   2.74
    29 0.200E+05 0.269E-03 0.425E+01 0.544E-04   1.00      0  -0.00   2.48
itn     obj       ng        delt      ns        rho    trit   phi    dn   
    30 0.200E+05 0.958E-04 0.425E+01 0.142E-04   0.47      0  -0.00   2.32

************************ CUTEst statistics ************************

 Code used               :  LDLTR
 Problem                 :  BDQRTIC   
 # variables             =            5000
 # objective functions   =          102.00
 # objective gradients   =          102.00
 # objective Hessians    =            0.00
 # iterations            =              30
 Exit code               =               0
 Final f                 =   2.0006257E+04
 Final ||g||             =   9.5808922E-05
 Set up time             =            0.05 seconds
 Solve time              =            7.72 seconds

******************************************************************
```


### Matlab 
To try-out an example in Matlab you can navigate to MATLAB/.
On the command prompt you can then issue the command

	>> EXAMPLE
	
The outcomes of running EXAMPLE.m look like:
```
>> EXAMPLE
 Solver                          :          LDLtr
 Line search                     :          Wolfe

  Itn    Objective    Norm g     Delta    Norm s    Rho      TRit   Phi     min(D)
    0 2.9115000e+03  1.53e+03  9.19e+00  4.59e+00    --    --     --        --    
    1-8.2765519e+02  4.73e+01  9.19e+00  5.35e+00  9.97e-01 0   0.00e+00 1.12e+02
    2-8.4272248e+02  2.29e+01  9.19e+00  4.70e-01  1.36e+00 0   0.00e+00 1.12e+02
    3-8.4812344e+02  6.35e+00  9.19e+00  4.34e-01  1.13e+00 0   0.00e+00 8.31e+01
    4-8.4841355e+02  2.47e+00  9.19e+00  8.81e-02  1.22e+00 0   0.00e+00 8.13e+01
    5-8.4847469e+02  1.43e+00  9.19e+00  3.82e-02  1.39e+00 0   0.00e+00 8.18e+01
    6-8.4849823e+02  5.31e-01  9.19e+00  3.15e-02  1.16e+00 0   0.00e+00 7.09e+01
    7-8.4849994e+02  1.07e-01  9.19e+00  7.01e-03  1.13e+00 0   0.00e+00 6.89e+01
    8-8.4850000e+02  6.98e-03  9.19e+00  1.17e-03  1.05e+00 0   0.00e+00 6.53e+01
    9-8.4850000e+02  5.52e-04  9.19e+00  1.23e-04  1.05e+00 0   0.00e+00 6.74e+01
   10-8.4850000e+02  8.23e-05  9.19e+00  8.55e-06  1.09e+00 0   0.00e+00 6.67e+01

 PROBLEM EXPERIM   --- Optimal                 : Itns  =    10, Num f =    23, Cpu time   = 0.007, n     =    5

 STATS    EXPERIM   : Time  = 0.007, Skipped =      0

```

