# LDLtr
                                                                       
 LDLtr: An LDLT trust-region region algorithm for the minimization 
 of a general nonlinear function when only gradients are available.                                     
                                                                       
                            minimize f(x)                                   
 
 The methods seek a point "norm(grad(x)) = 0", where grad(x) is
 the gradient of the objective f(x).
 
## Setup
Files are organized in the following:

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
from FORTRAN/cutest/src/ldltr to "LOCAL_CUTEST"/cutest/scr/.
When LDLtr is interfaced with CUTEst one can run the algorithm with e.g.,

		runcutest -p ldlqn -D BDQRTIC

where BDQRTIC is a problem name from CUTEst. 

### Matlab
With CUTEst for Matlab installed one can modify the path specified
in MATLAB/runUCset.m to point ot the local CUTEst installation.

## Examples

### Fortran
To try-out an example 


### Matlab 
To try-out three examples navigate inside folder "EXPERIMENTS/".
From within the folder you can run the examples by

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

