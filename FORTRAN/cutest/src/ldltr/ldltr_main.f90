PROGRAM LDLTR_main

!  LDLTR test driver for problems derived from SIF files.

!  Johannes J. Brust, Nov. 2023

      use ldlqn

      INTEGER :: n, derivs, maxit, i, status, statusLDLTR, ibound, lwa, nf, ng, nh, itn
      INTEGER :: idiff, kmax
      INTEGER, PARAMETER :: input = 55, out = 6, inspec = 56
      INTEGER, PARAMETER :: io_buffer = 11
      DOUBLE PRECISION :: biginf, eps, fopt, gnorm
      DOUBLE PRECISION, PARAMETER :: tinyL = 1.0D-6
      CHARACTER ( LEN = 10 ) :: pname
      DOUBLE PRECISION :: CPU( 2 ), CALLS( 4 )
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: JBOUND
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: X, BL, BU, WA, g
      CHARACTER ( LEN = 10 ), ALLOCATABLE, DIMENSION( : )  :: XNAMES
      EXTERNAL :: feval, geval

!  Open the spec file

      OPEN( inspec, FILE = 'LDLTR.SPC', FORM = 'FORMATTED',                    &
            STATUS = 'OLD' )
      REWIND( inspec )

!  read input Spec data

! DERIVS      <--  Derivatives available (<=0 none, 1 1st, >=2 1st and 2nd)
! BIGINF      <--  Bounds larger than biginf in magnitude are infinite
! EPS         <--  Stopping tolerance
! NF          <--  Maximum number of function evaluations (0 -> 1000*n)
! IDIFF       <--  Forward (<=1) or central (>=2) differences
! KMAX        <--  Hessian recomputation interval

!     READ( inspec, "( I10 )" ) derivs
!     READ( inspec, "( E10.3 )" ) biginf
!     READ( inspec, "( E10.3 )" ) eps
!     READ( inspec, "( I10 )" ) nf
!     READ( inspec, "( I10 )" ) idiff
!     READ( inspec, "( I10 )" ) kmax

!       READ( inspec, "( I10, /, E10.3, /, E10.3, /, I10, /, I10, /, I10 )" ) &
!         derivs, biginf, eps, nf, idiff, kmax
! !     READ( inspec, "( I10, /, 2( E10.3, /), (2I10, / ), I10 )" )              &
! !       derivs, biginf, eps, nf, idiff, kmax
!       idiff = MIN( MAX( 1, idiff ), 2 )
!       kmax = MAX( kmax, 1 )

!  close input file

      CLOSE( inspec )

!  open the input data file

      OPEN( input, FILE = 'OUTSDIF.d', FORM = 'FORMATTED',                     &
            STATUS = 'OLD' )
      REWIND( input )

!  find the problem dimension

      CALL CUTEST_udimen( status, input, n )
      IF ( status /= 0 ) GO TO 910

!  allocate workspace
      
      ALLOCATE( X( n ), g( n ), BL(n), BU( n ), XNAMES( n ),              &
                STAT = status )
      IF ( status /= 0 ) GO TO 990

!  set up SIF data

      CALL CUTEST_usetup( status, input, out, io_buffer, n, X, BL, BU )

!  obtain variable names

      CALL CUTEST_unames( status, n, pname, XNAMES )
      IF ( status /= 0 ) GO TO 910

!  call the optimizer

  call ldlqn_bfgswtr(X, itn, nf, statusLDLTR)

!  output solution

  CALL CUTEST_ureport( status, CALLS, CPU )
  IF ( status /= 0 ) GO TO 910

! call the objective function and gradient at the computed point

  CALL feval(n,X,fopt)
  CALL geval(n,X,g)

      
  gnorm = norm2(g)
      
      WRITE ( out, "( /, 24('*'), ' CUTEst statistics ', 24('*') //,           &
     &         ' Code used               :  LDLTR',   /,                       &
     &         ' Problem                 :  ', A10,    /,                      &
     &         ' # variables             =      ', I10, /,                     &
     &         ' # objective functions   =        ', F8.2, /,                  &
     &         ' # objective gradients   =        ', F8.2, /,                  &
     &         ' # objective Hessians    =        ', F8.2, /,                  &
     &         ' # iterations            =      ', I10, /,                     &
     &         ' Exit code               =      ', I10 /,                      &
     &         ' Final f                 = ', ES15.7, /,                       &
     &         ' Final ||g||             = ', ES15.7, /,                       &
     &         ' Set up time             =      ', 0P, F10.2, ' seconds', /,   &
     &         ' Solve time              =      ', 0P, F10.2, ' seconds', //,  &
     &           66('*') / )" ) pname, n, ( CALLS( i ), i = 1, 3 ),            &
                   itn, statusLDLTR, fopt, gnorm, CPU( 1 ), CPU( 2 )
      CLOSE( input  )
      STOP

  910 CONTINUE
      WRITE( out, "( ' CUTEst error, status = ', I0, ', stopping' )" ) status
      STOP

  990 CONTINUE
      WRITE( out, "( ' Allocation error, status = ', I0 )" ) status
      STOP
      END

      SUBROUTINE feval( n, X, f )
      INTEGER :: n
      DOUBLE PRECISION :: f, X( n )
      INTEGER :: status
      INTEGER, PARAMETER :: out = 6
      CALL CUTEST_ufn( status, n, X, f )
      IF ( status /= 0 ) THEN
        WRITE( out, "( ' CUTEst error, status = ', i0, ', stopping' )" )       &
          status
        STOP
      END IF
      RETURN
      END

      SUBROUTINE geval( n, X, G )
      INTEGER :: n
      DOUBLE PRECISION :: X( n ), G( n )
      INTEGER :: status
      INTEGER, PARAMETER :: out = 6
      CALL CUTEST_ugr( status, n, X, G )
      IF ( status /= 0 ) THEN
        WRITE( out, "( ' CUTEst error, status = ', i0, ', stopping' )" )       &
          status
        STOP
      END IF
      RETURN
      END

      
