
!  Dummy LDLTR for testing LDLTR_main interface to CUTEst

   SUBROUTINE ldlqn_bfgswtr(X0, itn, nf, status)
   integer :: itn, nf, status   
   double precision :: X0( : )
   itn      = 0
   nf       = 0
   status   = 0 
   X0(1)    = 0.0
   END SUBROUTINE ldlqn_bfgswtr

   