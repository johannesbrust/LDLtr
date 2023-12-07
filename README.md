# LDLtr
                                                                       
 LDLtr: An LDLT trust-region region algorithm for the minimization 
 of a general nonlinear function when only gradients are available.                                     
                                                                       
                            minimize f(x)                                   
 
 The methods seek a point "norm(grad(x)) = 0", where grad(x) is
 the gradient of the objective f(x).
 
## Setup
Files are organized in 6 folders:

	ALGS/ 		(main algorithms)
	AUXILIARY/ 	(support routines)
	DATA/ 		(data for figures and tables)
	EXPERIMENTS/ 	(examples and experiments)
	EXTERNAL/ 	(third party software)
	FIGS/		(figures)

To rerun any experiment from the folder "EX_COMP_LMSSM_LSR1_m_q5" the 
CUTEst (https://github.com/ralna/CUTEst/wiki) test set has to be installed. 
If CUTEst is installed, then the relative path to the CUTEst installation
is to be updated in 

	AUXILIARY/CUTEst_init.m

Without a CUTEst installation, the figures can still be plotted 
from within the "AUXILIARY/" folder by calling 
fig*.m, *={1,2,3,4,5,6,7,8}

Detailed descriptions of the signatures of the algorithms in "ALGS/"
are found in the beginning comments of the respective files

The codes have been tested on:

Matlab R2016a, macOS Catalina, {CUTEst n/a, July 2017}

## Examples

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

