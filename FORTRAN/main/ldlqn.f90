!>
!> ldlqn.f90
!> 
!> Implementation of the Matlab algorithm bfgsWTRL2SS9.m
!> This method updates a LDLT factorization of a quasi-Newton matrix
!> and uses a trust-region strategy to compute search directions
!>
!> 
!>-------------------------------------------------------------------
!> 09/11/23, J.B., initial implementation

module ldlqn

    use kind_parameter, only : dp
    use ldlqn_bfgsRDRup
    use ldlqn_invLDLCG
    use ldlqn_triang
    use ldlqn_pars
    
    implicit none

    contains

    subroutine ldlqn_bfgswtr(x, itn, nf, status) ! , feval, geval, 
        !> Comments of original Matlab implementation
        ! %bfgsWTRL2SS9  [x,itn,nf,skipped,status] = bfgsWTRL2SS9(prob,outfile) finds
        ! %        a local minimizer of a scalar-valued multivariate function
        ! %        using a variant of the BFGS quasi-Newton method with a
        ! %        weighted trust-region subproblem.
        ! %
        ! %        prob is a structure that provides the details of the
        ! %        CUTEst problem being solved.
        ! %
        ! %        On successful termination, x is the first point found at
        ! %        which the norm of the gradient is less than a tolerance.
        ! %
        ! %        Input arguments are set in the file bfgs_parms.m.
        ! %
        ! %        The meaning of each line of output is as follows
        ! %
        ! %        Itn          The current iteration number.
        ! %
        ! %        f(x)         The current value of the objective function f(x).
        ! %
        ! %        Norm g       The current value of the gradient of f(x). In
        ! %                     the neighborhood of a local minimizer, Norm g
        ! %                     will be small.
        ! %
        ! %        Delta        The radius of the trust-region.
        ! %
        ! %        Rho          Sufficient decrease metric: (fk1-fk)/Q(sk)
        ! %
        ! %        TRit         Trust-region subproblem iterations
        ! %
        ! %        Phi          Subproblem constraint: |1/Delta - 1/norm(sk)|
        ! %
        ! %--------------------------------------------------------------------------
        ! % This version extends an implementation from 29 Feb 2020 by Philip E. Gill
        ! %
        ! % The new version is first implemented 24 Jan 2023 by Johannes J. Brust
        ! % University of California, San Diego.
        ! %
        ! % Updates: 
        ! % 01/26/23, J.B., Initial version of WTRL2
        ! % 02/10/23, J.B., Scaled trust-region subproblem
        ! % 02/14/23, J.B., Updated scaled problem
        ! % 02/16/23, J.B., Remove normalization
        ! % 02/23/23, J.B., Add diagonal shift in subproblem
        ! % 02/24/23, J.B., Low rank approximation in subproblem
        ! % 03/01/23, J.B., Removing function evaluations
        ! % 03/02/23, J.B., Use of inverse update, store and update only inverse
        ! % 05/09/23, J.B., Function implements a CG solver in the subproblem
        ! % 05/16/23, J.B., Further correction of CG in Newton solver
        ! % 05/23/23, J.B., Use of a banded approximation for the TR subproblem


        !>
        !> Fortran implementation comments
        !>    
        !> The input
        !>    
        !>  feval(n,x,f)    -- Objective function
        !> 
        !> needs to be defined in the main program. n is the problem size, x is 
        !> current solution estimate and f is computed function value.
        !>    
        !>  geval(n,x,g)    -- Objective gradient
        !> 
        !> is the function that evaluates the gradient at x. This function
        !> also needs to be implemented in the main program.
        !> 
        !>  x               -- Solution estimate (size(x)=n)
        !>
        !>  itn             -- Number of iterations
        !>
        !>  nf              -- Number of function evaluations
        !>
        !>  status          -- Flag for the outcome of a run 
        !>

        !>
        !> Initializations
        !>

        real(dp), intent(out)   ::          x(:)
        integer, intent(out)    ::          itn, nf, status 

        integer                 ::          n, solveTRSUB, trIt, i, hasGRADMIN, &
                                            itMin, numAccept, numTRInc, numTRDec, &
                                            skipped, optimal, unbounded, nearoptimal, &
                                            itnlimit, badsearch, info

        real(dp)                ::          rho, phi, Mjj, MjjMax, MjjMin, &
                                            normg, normgMax, fMax, dn, f, f1, f1t, f2t, &
                                            nsk, nyk, sy, delta, sig, ared, &
                                            ng1, relphi, sigCG, fCG, ngcg, rho1, &
                                            pred, pred1, ared1, ng1t, skskp, ng2t, relphiN, &
                                            sigN, sigb, phib, siga 

        real(dp), allocatable   ::          ons(:), D(:), DI(:), R(:,:), GI(:), &
                                            g(:), g1(:), g1t(:), g2t(:), sk(:), sk1(:), yk(:), b1(:), b2(:), &
                                            vk(:), vkp(:), skp(:), hk(:), gCG(:), x1t(:), &
                                            x1(:), skN(:)

        real(dp), allocatable   ::          dbgv1(:), dbgv2(:)

        n                   = size(x)
        solveTRSUB          = 0
        trIt                = 0
        i                   = 1
        hasGRADMIN          = 0
        itMin               = 1
        numAccept           = 0
        numTRInc            = 0
        numTRDec            = 0
        skipped             = 0
        optimal             = 0
        unbounded           = 0
        nearoptimal         = 0
        itnlimit            = 0
        badsearch           = 0
        status              = 0
        rho                 = 0.0
        phi                 = 0.0
        
        MjjMax              = 1e4
        MjjMin              = 1e-2

        allocate(ons(n),    source=1.0_dp)
        allocate(D(n),      source=1.0_dp)
        allocate(DI(n),     source=1.0_dp)
        allocate(R(n,n),    source=0.0_dp)
        allocate(GI(n),     source=1.0_dp)
        allocate(g(n),      source=1.0_dp)
        allocate(g1(n),     source=1.0_dp)
        allocate(g1t(n),    source=1.0_dp)
        allocate(g2t(n),    source=1.0_dp)
        allocate(sk(n),     source=0.0_dp)
        allocate(sk1(n),    source=0.0_dp)
        allocate(skN(n),    source=0.0_dp)
        allocate(yk(n),     source=0.0_dp)
        allocate(b1(n),     source=0.0_dp)
        allocate(b2(n),     source=0.0_dp)
        allocate(vk(n),     source=0.0_dp)
        allocate(vkp(n),    source=0.0_dp)
        allocate(skp(n),    source=0.0_dp)
        allocate(hk(n),     source=0.0_dp)
        allocate(gCG(n),    source=1.0_dp)
        allocate(x1t(n),    source=0.0_dp)
        allocate(x1(n),     source=0.0_dp)
        

        !
        ! For debugging
        !
        ! allocate(dbgv1(n),     source=0.0_dp)
        ! allocate(dbgv2(n),     source=0.0_dp)

        !> 
        !> Initial call to objective and gradient
        !>
        call feval(n,x,f)
        call geval(n,x,g)
        nf          = nf + 1


        !>
        !> Initializations of the quasi-Newton factorization
        !> 
        normg       = norm2(g)
        normgMax    = 1 + normg 
        fMax        = 1 + abs(f)

        Mjj         = 1/max(1.0,normg)
        Mjj         = min(max(Mjj,MjjMin),MjjMax)
        DI          = Mjj * DI
        D           = 1.0 / DI
        dn          = minval(D)

        do i = 1,n 
            R(i,i) = 1.0
        end do

        if (print_run_det == 1) then
            print *, '************************************************************************'
            print *, '*'
            print *, '* LDLTR1'
            print *, '* '
            print *, '* LDL Trust-Region Quasi-Newton Algorithm'
            print *, '*'
            print ' (a5,i6)', '* n=', n
            print ' (a6,f6.3)',' * fx=', f
            print ' (a6,f6.3)',' * gx=', normg
            print *, '*'
            print *, '************************************************************************'
            print *, ''            
            !print '(a6,a6,a6,a6,a6,a6,a6,a6,a6)', 'itn   ', ' obj  ', ' ng   ', ' delt ', &
            !        ' ns   ', ' rho  ', ' trit ', ' phi  ', ' dn   '
            ! print '(a6,a9,a9,a9,a9,a9,a6,a9,a9)', 'itn   ', ' obj     ', ' ng      ', ' delt    ', &
            !         ' ns      ', ' rho     ', ' trit ', ' phi     ', ' dn      '
            print '(a6,1x,a9,1x,a9,1x,a9,1x,a9,1x,a6,1x,a6,1x,a6,1x,a6)', 'itn   ', ' obj     ', ' ng      ', ' delt    ', &
                    ' ns      ', ' rho  ', ' trit ', ' phi  ', ' dn   '                    
        end if

        !> 
        !> Use of an initial line-search to start the method
        !> 

        x1(:)   = x(:)
        x1t(:)  = x(:)        
        g1(:)   = g(:)
        b1(:)   = g(:)
        sk      = -Mjj * g
        nsk     = norm2(sk)
        f1      = f

        stpmax  = dxInf/ (1+nsk)
        stp     = dxMax/ (1+nsk)
        stp     = min(1.0,stp)
        stp     = min(stpmax,stp)

        call cvsrch(n,x1,f1,g1,sk,stp,ftol,gtol,xtol,stpmin,stpmax,maxfev,info,nf,b1)
        
        !
        ! For debugging. Print 1st step
        !
        ! print '(a3,i6)','nf:', nf
        ! print '(a4,e10.3)','stp:', stp
        ! print '(a3,e10.3)','ng:', norm2(g1)
        ! print '(a2,e10.3)','f:', f1

        if (info > 0 .and. info < 6) then

            sk      = x1 - x
            yk      = g1 - g
            nsk     = norm2(sk)
            nyk     = norm2(yk)
            sy      = dot_product((sk/nsk),yk/nyk) * nsk * nyk

            ! x       = x1
            ! f       = f1
            ! g       = g1
            ! normg   = norm2(g)

        end if

        !>
        !> Quasi-Newton update
        !>

        if (sy > 0) then

            !> Product Hk*yk, where Hk is the factored inverse 
            !> Quasi-Newton matrix. The vector b2 stores the result

            call lowermult(transpose(R),yk,b1)
            b1 = DI*b1
            call uppermult(R,b1,b2)

            !> Two rank-1 updates to the factorization
            call bfgsRDRupEig(R, DI, GI, sk, yk, b2, nsk, nyk ) 

            D   = 1.0 / DI
            dn  = minval(D)

        end if 

        ! Setting the initial trust-region radius
        delta   = 2 * nsk

        !>
        !> Main loop
        !>

        do itn = 1, maxIterations

            !> Check for status
            if (normg <= tolStat) then
                optimal = 1
            else if (f <= unBoundedf) then                
                unbounded = 1
            else if (      (count_near == 1 .and. badsearch == 1) .and. &
                            ((normg <= (2e-16**(2.0/3.0))*normgMax) .or. & 
                            abs(f) <= (2e-16**(2.0/3.0))*fMax)          &
                    ) then
                nearoptimal = 1     
            end if

            !>
            !> Printing results
            !>

            if (print_run_det == 1) then

                if (mod(itn,lines) == 0) then
                    print '(a6,1x,a9,1x,a9,1x,a9,1x,a9,1x,a6,1x,a6,1x,a6,1x,a6)', 'itn   ', ' obj     ', ' ng      ', ' delt    ', &
                    ' ns      ', ' rho  ', ' trit ', ' phi  ', ' dn   '

                    ! print '(a6,a9,a9,a9,a9,a9,a6,a9,a9)', 'itn   ', '  obj    ', '  ng     ', '  delt   ', &
                    ! '  ns     ', '  rho    ', '  trit', '  phi    ', '  dn     '
                end if

                ! rho     = 0.0
                ! trIt    = 0
                ! phi     = 0.0
                
                print '(i6,1x,e9.3,1x,e9.3,1x,e9.3,1x,e9.3,1x,f6.2,1x,i6,1x,f6.2,1x,f6.2)', itn, f, normg, delta, nsk, &
                        rho, trIt, phi, dn

            end if

            if (optimal == 1 .or. unbounded == 1 .or. nearoptimal == 1 .or. & 
                itn == maxIterations .or. badsearch == 1) then
                exit
            endif

            !>
            !> Computing the unconstrained minimizer
            !> 

            call lowermult(transpose(R),g,hk)

            sig = 0.0
            vk  = -( hk / (D+sig) )
            b1  = -( DI * hk )
            
            call uppermult(R,b1,sk)        ! Unconstrained minimizer sk

            nsk = norm2(sk)

            !
            ! For debugging. Print norm of unconstrained step
            !            
            ! print '(a4,e10.3)','nsk:', nsk
            ! print '(a5,e10.3)','Delk:', delta


            if (nsk - delta > 0) then      ! Check if solving the subproblem is needed
                solveTRSUB = 1
            else
                solveTRSUB = 0
            end if

            sk1(:)  = sk(:)                   ! Store the current best step
            b2(:)   = sk(:)                   ! And store the quasi-Newton step

            call feval(n,x+sk,f1t)         ! Evaluate function and gradient at step 
            call geval(n,x+sk,g1t)
            nf      = nf + 1
            ared1   = f - f1t
            ng1t    = norm2(g1t)

            !>
            !> In case quasi-Newton step is optimal terminate early
            !>
            if (ng1t <= tolStat .and. delta < tolStat) then 
                solveTRSUB = 0
            end if

            trIt = 0
            
            !>
            !> Phase 1 of algorithm. Solving a diagonal trust-region
            !> subproblem. 
            !>
            if (solveTRSUB == 1) then       ! .and. nmax < n
                                            ! Branching strategy is not implemented in this version

                phi     = 1/nsk - 1/delta
                relphi  = (delta-nsk)/ delta

                vkp     = -(GI * vk) / (D + sig* GI)

                call uppermult(R,vkp,skp)

                !> 
                !> Trust-region loop
                !> 

                do i = 1, trMax 

                    ! Terminate the subproblem if a "solution" on the boundary is found
                    if (abs(relphi) <= trTol) then        
                        exit
                    end if

                    !>
                    !> 1D Newton's iteration
                    !>

                    skskp   = dot_product(sk,skp)
                    sig     = sig + ( (nsk*nsk) / skskp ) * ((delta-nsk)/ delta)

                    !vk      = -(GI * hk) / (D + sig* GI)
                    vk      = -(hk) / (D + sig* GI)

                    call uppermult(R,vk,sk)

                    ! Further checks on function improvement
                    ! Performed when the gradient is already relatively small
                    if (normg < 1) then
                        
                        call feval(n,x+sk,f2t)
                        call geval(n,x+sk,g2t)
                        nf = nf + 1
                        
                        if (f-f2t > ared1 .and. sig > 0) then 

                            ared1   = f-f2t
                            sk1     = sk
                            f1t     = f2t
                            g1t     = g2t

                        end if

                        if (norm2(g2t) < tolStat .and. delta < tolStat) then 

                            ared1   = f-f2t
                            sk1     = sk
                            f1t     = f2t
                            g1t     = g2t

                            exit

                        end if

                    end if

                    nsk     = norm2(sk)

                    phi     = 1/nsk - 1/delta
                    relphi  = (delta-nsk)/ delta

                    vkp     = -(GI * vk) / (D + sig* GI)

                    call uppermult(R,vkp,skp)

                    trIt    = trIt + 1

                end do
                
                ! Bi-section method is called if Newton's method failed

                skN     = sk
                relphiN = relphi
                sigN    = sig

                if (sig < 0 .or. trIt == trMax) then

                    !>
                    !> Bisection algorithm
                    !>

                    sk      = b2            ! Restore the quasi-Newton step         
                    nsk     = norm2(sk)

                    ! For the next parts vk -> vkb and vkp -> skb 
                    ! when compared to the Matlab's implementation naming convention

                    sigb    = 10*maxval(D)             ! Taking a reasonably large sigma value
                    vk      = - (hk) / (D + sigb*DI)
                    
                    call uppermult(R,vk,vkp)

                    phib    = 1/norm2(vkp) - 1/delta

                    ! Search for upper sigma estimate to bracket the 
                    ! secular equation (required for the bisection method) 
                    do i = 2,10

                        if (phib > 0) then 
                            exit
                        end if
                        
                        sigb    = (10**i) * maxval(D)
                        vk      = - (hk) / (D + sigb*DI)
                    
                        call uppermult(R,vk,vkp)

                        phib    = 1/norm2(vkp) - 1/delta

                    end do

                    ! Preparation for the Bi-section iteration

                    siga    = 0.0
                    phi     = 1/nsk - 1/delta
                    relphi  = (delta-nsk) / delta

                    !>
                    !> Bisection iteration
                    !> 
                    do i = trIt, 2*trMax

                        ! Terminate early
                        if (abs(relphi) <= trTol ) then 
                            exit
                        end if

                        sig     = (siga + sigb) / 2 

                        vk      = - (hk) / (D + sigb*DI)
                    
                        call uppermult(R,vk,sk)

                        nsk     = norm2(sk)

                        phi     = 1/nsk - 1/delta

                        relphi  = (delta-nsk) / delta


                        ! Further checks on function improvement
                        ! Performed when the gradient is already relatively small
                        ! Analogous to the Newton iteration checks
                        if (normg < 1) then
                            
                            call feval(n,x+sk,f2t)
                            call geval(n,x+sk,g2t)
                            nf = nf + 1
                            
                            if (f-f2t > ared1 .and. sig > 0) then 

                                ared1   = f-f2t
                                sk1     = sk
                                f1t     = f2t
                                g1t     = g2t

                            end if

                            if (norm2(g2t) < tolStat .and. delta < tolStat) then 

                                ared1   = f-f2t
                                sk1     = sk
                                f1t     = f2t
                                g1t     = g2t

                                exit

                            end if

                        end if

                        if (phi < 0) then 
                            siga = sig
                        else 
                            sigb = sig 
                        end if 
                    
                    end do

                end if ! End of the bisection part

                ! Pick whichever root solver performed better
                if (sigN >= 0 .and. abs(relphiN) <= abs(relphi) ) then 

                    sk  = skN
                    sig = sigN

                end if

                ng1t = norm2(g1t)

            end if ! End of trust-region subproblem

            !>
            !> Phase 2 of the Algorithm. This computes a sequence of shift values combined
            !> with a conjugate gradient iteration
            !>

            sigCG       = sig
            call feval(n,x+sk,fCG)
            call geval(n,x+sk,gCG)
            nf          = nf + 1
            ngcg        = norm2(gCG)

            itMin       = 0
            hasGRADMIN  = 0

            if (solveTRSUB == 1) then 

                do i = 1, maxcg

                    call invLDLCG(vk, b1, b2, sigCG, -hk, R, D, tolICG, maxitICG )

                    call uppermult(R,vk,skN)

                    call feval(n,x+skN,f2t)
                    call geval(n,x+skN,g2t)
                    nf          = nf + 1
                    ng2t        = norm2(g2t)

                    ! Selecting a step from the CG iteration
                    if (f2t < fCG .and. abs(f2t-f) / abs(f) > 1e-12) then
                        
                        sk      = skN
                        fCG     = f2t
                        itMin   = i 

                    else if (abs(f2t-f)/abs(f) < 1e-12 .and. ng2t < ngcg & 
                                    .and. delta < 1e-7) then 

                        sk          =  skN
                        fCG         = f2t 
                        ngcg        = ng2t 
                        itMin       = i
                        hasGRADMIN  = 1

                    end if

                    ! Shrinking the shift parameter
                    sigCG   = shrink * sigCG

                end do

            end if

            ! Comparing ratios
            ! First the ratio based on the best step from phase 1 is computed
            ! Second the ratio based on the step from phase 2 is computed (before comparing them)
            call backsub(R,sk1,b1)
            
            call forwardsub(transpose(R),D*b1,b2)

            pred1   = -(dot_product(g,sk1) + 0.5 * dot_product(b2,sk1) )

            rho1    = ared1 / pred1

            ! Early termination
            if (ng1t <= tolStat .and. delta < tolStat) then 
                rho1 = 1
            end if

            ! If model is not positive definite, omit the step
            if (pred1 < 0.0) then 
                sk1(:)  = 0.0
                rho1    = -1.0
            end if

            !>
            !> Computing a trial step
            !> And evaluating it
            !> 

            x1t         = x + sk
            call feval(n,x1t,f2t)
            call geval(n,x1t,g2t)
            nf          = nf + 1

            ared        = f - f2t
          
            call backsub(R,sk,b1)
            
            call forwardsub(transpose(R),D*b1,b2)

            pred    = -(dot_product(g,sk) + 0.5 * dot_product(b2,sk) )

            rho     = ared / pred

            !
            ! For debugging. Print rho1 and rho
            !
            ! print '(a4,e10.3)','nsk1:', norm2(sk1)
            ! print '(a4,e10.3)','nsk:', norm2(sk)            
            ! print '(a4,e10.3)','rho1:', rho1
            ! print '(a5,e10.3)','rho:', rho
            ! print '(a5,e10.3)','ared1:', ared1
            ! print '(a5,e10.3)','pred1:', pred1
            ! print '(a5,e10.3)','ared:', ared
            ! print '(a5,e10.3)','pred:', pred
            ! print '(a5,e10.3)','sig:', sig
            ! print '(a5,i6)','itMin:', itMin
            ! print '(a5,e10.3)','nx1t:', norm2(x1t)
            ! print '(a5,e10.3)','f:', f
            ! print '(a5,e10.3)','f2t:', f2t
            ! print '(a5,e10.3)','sum(x):', sum(x)
            ! print '(a5,i6)','i:', i

            ! If model is not positive definite, omit the step
            if (pred < 0.0) then 
                sk(:)  = 0.0
                rho    = -1.0
            end if

            ! Update the relevant quantities based on the current 
            ! best point
            if (rho1 > rho) then
                
                rho     = rho1
                sk      = sk1 
                x1t     = x + sk
                
            else 

                f1t     = f2t
                g1t     = g2t 

            end if

            yk      = g1t - g
            nsk     = norm2(sk)
            nyk     = norm2(yk)

            !>
            !> Trust-region updating
            !> 

            ! Accepting a step
            if (c1 < rho) then 

                x1  = x1t
                g1  = g1t
                f1  = f1t
                
                numAccept = numAccept + 1

                ! Shrinking value
                if (itMin == maxcg) then
                    ! Halve the shrink value 
                    shrink = 0.5*shrink
                elseif (itMin <= 1) then 
                    ! Double the shrink value
                    shrink = (1.0/0.5)*shrink
                end if
                shrink = min(shrink,0.25)
                shrink = max(shrink,0.25**10)

            else 
                ! Step is not accepted
                x1  = x
                g1  = g
                f1  = f

            end if

            ! Updating the radius
            if (c2 < rho) then 

                if (nsk <= c3*delta) then 
                    delta = 1.0*delta
                else 
                    delta = c4*delta
                    numTRInc = numTRInc + 1
                end if

            elseif (c5 <= rho .and. rho <= c6) then 
                delta = 1.0*delta
            else
                delta = c7*delta
                numTRDec = numTRDec + 1
            end if 

            !>
            !> BFGS update
            !>

            sy          = dot_product(sk,yk)
            normg       = norm2(g1)
            normgMax    = max(normgMax,normg) 

            if (sy > 0 .and. rho >= 0) then 

                !> Product Hk*yk, where Hk is the factored inverse 
                !> Quasi-Newton matrix. The vector b2 stores the result

                call lowermult(transpose(R),yk,b1)
                b1 = DI*b1
                call uppermult(R,b1,b2)

                !> Two rank-1 updates to the factorization
                call bfgsRDRupEig(R, DI, GI, sk, yk, b2, nsk, nyk ) 

                D   = 1.0 / DI
                dn  = minval(D)

                !
                ! Debug: Secant condition
                !
                ! call lowermult(transpose(R),yk,dbgv1)
                ! dbgv1 = DI*dbgv1
                ! call uppermult(R,dbgv1,dbgv2)

                ! print '(a5,e10.3)','Sec:', norm2(dbgv2-sk)



            end if 

            x = x1
            f = f1
            g = g1

            if (delta < 1e-22) then                
                badsearch = 1
            end if

        end do

        if (optimal == 1) then 
            status = 0
        elseif (unbounded == 1) then 
            status = 1
        elseif (nearoptimal == 1) then 
            status = 2
        elseif (itn == maxIterations) then
            status = 3
        elseif (badsearch == 1) then
            status = 4
        endif


    end subroutine

end module