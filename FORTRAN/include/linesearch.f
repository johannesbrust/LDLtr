c
c     ----------------------------------------------------------------------
c     Line Search Routine
c     ----------------------------------------------------------------------
c
      subroutine cvsrch(n,x,f,g,s,stp,ftol,gtol,xtol,
     *                  stpmin,stpmax,maxfev,info,nfev,wa)
      integer n,maxfev,info,nfev
      double precision f,stp,ftol,gtol,xtol,stpmin,stpmax
      double precision x(n),g(n),s(n),wa(n)
c     **********
c
c     Subroutine cvsrch
c
c     The purpose of cvsrch is to find a step which satisfies 
c     a sufficient decrease condition and a curvature condition.
c     The user must provide a subroutine which calculates the
c     function and the gradient.
c
c     At each stage the subroutine updates an interval of
c     uncertainty with endpoints stx and sty. The interval of
c     uncertainty is initially chosen so that it contains a 
c     minimizer of the modified function
c
c          f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).
c
c     If a step is obtained for which the modified function 
c     has a nonpositive function value and nonnegative derivative, 
c     then the interval of uncertainty is chosen so that it 
c     contains a minimizer of f(x+stp*s).
c
c     The algorithm is designed to find a step which satisfies 
c     the sufficient decrease condition 
c
c           f(x+stp*s) .le. f(x) + ftol*stp*(gradf(x)'s),
c
c     and the curvature condition
c
c           abs(gradf(x+stp*s)'s)) .le. gtol*abs(gradf(x)'s).
c
c     If ftol is less than gtol and if, for example, the function
c     is bounded below, then there is always a step which satisfies
c     both conditions. If no step can be found which satisfies both
c     conditions, then the algorithm usually stops when rounding
c     errors prevent further progress. In this case stp only 
c     satisfies the sufficient decrease condition.
c
c     The subroutine statement is
c
c        subroutine cvsrch(fcn,n,x,f,g,s,stp,ftol,gtol,xtol,
c                          stpmin,stpmax,maxfev,info,nfev,wa)
c     where
c
c	fcn is the name of the user-supplied subroutine which
c         calculates the function and the gradient.  fcn must 
c      	  be declared in an external statement in the user 
c         calling program, and should be written as follows.
c
c	  subroutine fcn(n,x,f,g)
c         integer n
c         f
c         x(n),g(n)
c	  ----------
c         Calculate the function at x and
c         return this value in the variable f.
c         Calculate the gradient at x and
c         return this vector in g.
c	  ----------
c	  return
c	  end
c
c       n is a positive integer input variable set to the number
c	  of variables.
c
c	x is an array of length n. On input it must contain the
c	  base point for the line search. On output it contains 
c         x + stp*s.
c
c	f is a variable. On input it must contain the value of f
c         at x. On output it contains the value of f at x + stp*s.
c
c	g is an array of length n. On input it must contain the
c         gradient of f at x. On output it contains the gradient
c         of f at x + stp*s.
c
c	s is an input array of length n which specifies the
c         search direction.
c
c	stp is a nonnegative variable. On input stp contains an
c         initial estimate of a satisfactory step. On output
c         stp contains the final estimate.
c
c       ftol and gtol are nonnegative input variables. Termination
c         occurs when the sufficient decrease condition and the
c         directional derivative condition are satisfied.
c
c	xtol is a nonnegative input variable. Termination occurs
c         when the relative width of the interval of uncertainty 
c	  is at most xtol.
c
c	stpmin and stpmax are nonnegative input variables which 
c	  specify lower and upper bounds for the step.
c
c	maxfev is a positive integer input variable. Termination
c         occurs when the number of calls to fcn is at least
c         maxfev by the end of an iteration.
c
c	info is an integer output variable set as follows:
c	  
c	  info = 0  Improper input parameters.
c
c	  info = 1  The sufficient decrease condition and the
c                   directional derivative condition hold.
c
c	  info = 2  Relative width of the interval of uncertainty
c		    is at most xtol.
c
c	  info = 3  Number of calls to fcn has reached maxfev.
c
c	  info = 4  The step is at the lower bound stpmin.
c
c	  info = 5  The step is at the upper bound stpmax.
c
c	  info = 6  Rounding errors prevent further progress.
c                   There may not be a step which satisfies the
c                   sufficient decrease and curvature conditions.
c                   Tolerances may be too small.
c
c       nfev is an integer output variable set to the number of
c         calls to fcn.
c
c	wa is a work array of length n.
c
c     Subprograms called
c
c	user-supplied......fcn
c
c	MINPACK-supplied...cstep
c
c	FORTRAN-supplied...abs,max,min
c	  
c     Argonne National Laboratory. MINPACK Project. June 1983
c     Jorge J. More', David J. Thuente
c
c     **********
      integer infoc,j
      logical brackt,stage1
      double precision dg,dgm,dginit,dgtest,dgx,dgxm,dgy,dgym,
     *       finit,ftest1,fm,fx,fxm,fy,fym,p5,p66,stx,sty,
     *       stmin,stmax,width,width1,xtrapf
      data p5,p66,xtrapf/0.5,0.66,4.0/
      info = 0
      infoc = 1
c
c     Check the input parameters for errors.
c
      if (n .le. 0 .or. stp .le. 0.0 .or. ftol .lt. 0.0 .or. 
     *    gtol .lt. 0.0 .or. xtol .lt. 0.0 .or. stpmin .lt. 0.0 
     *    .or. stpmax .lt. stpmin .or. maxfev .le. 0) return
c
c     Compute the initial gradient in the search direction
c     and check that s is a descent direction.
c
      dginit = 0.0
      do 10 j = 1, n
         dginit = dginit + g(j)*s(j)
   10    continue
      if (dginit .ge. 0.0) then
         return
      endif
c
c     Initialize local variables.
c
      brackt = .false.
      stage1 = .true.
      nfev = 0
      finit = f
      dgtest = ftol*dginit
      width = stpmax - stpmin
      width1 = width/0.5d0
      do 20 j = 1, n
         wa(j) = x(j)
   20    continue
c
c     The variables stx, fx, dgx contain the values of the step, 
c     function, and directional derivative at the best step.
c     The variables sty, fy, dgy contain the value of the step,
c     function, and derivative at the other endpoint of
c     the interval of uncertainty.
c     The variables stp, f, dg contain the values of the step,
c     function, and derivative at the current step.
c
      stx = 0.0
      fx = finit
      dgx = dginit
      sty = 0.0
      fy = finit
      dgy = dginit
c
c     Start of iteration.
c
   30 continue
c
c        Set the minimum and maximum steps to correspond
c        to the present interval of uncertainty.
c
         if (brackt) then
            stmin = min(stx,sty)
            stmax = max(stx,sty)
         else
            stmin = stx
            stmax = stp + xtrapf*(stp - stx)
         endif
c
c        Force the step to be within the bounds stpmax and stpmin.
c
         stp = max(stp,stpmin)
         stp = min(stp,stpmax)
c
c        If an unusual termination is to occur then let 
c        stp be the lowest point obtained so far.
c
         if ((brackt .and. (stp .le. stmin .or. stp .ge. stmax))
     *      .or. nfev .ge. maxfev-1 .or. infoc .eq. 0
     *      .or. (brackt .and. stmax-stmin .le. xtol*stmax)) stp = stx
c
c        Evaluate the function and gradient at stp
c        and compute the directional derivative.
c
         do 40 j = 1, n
            x(j) = wa(j) + stp*s(j)           
 40      continue

c         For debugging purposes only 
c         print *, x

         call feval(N,x,f)
         call geval(N,x,g)
         nfev = nfev + 1
         dg = 0.0
         do 50 j = 1, n
            dg = dg + g(j)*s(j)
   50       continue
         ftest1 = finit + stp*dgtest
c
c        Test for convergence.
c
         if ((brackt .and. (stp .le. stmin .or. stp .ge. stmax))
     *      .or. infoc .eq. 0) info = 6
         if (stp .eq. stpmax .and. 
     *       f .le. ftest1 .and. dg .le. dgtest) info = 5
         if (stp .eq. stpmin .and. 
     *       (f .gt. ftest1 .or. dg .ge. dgtest)) info = 4
         if (nfev .ge. maxfev) info = 3
         if (brackt .and. stmax-stmin .le. xtol*stmax) info = 2
         if (f .le. ftest1 .and. abs(dg) .le. gtol*(-dginit)) info = 1
c
c        Check for termination.
c
         if (info .ne. 0) return
c
c        In the first stage we seek a step for which the modified
c        function has a nonpositive value and nonnegative derivative.
c
         if (stage1 .and. f .le. ftest1 .and. 
     *       dg .ge. min(ftol,gtol)*dginit) stage1 = .false.
c
c        A modified function is used to predict the step only if
c        we have not obtained a step for which the modified
c        function has a nonpositive function value and nonnegative 
c        derivative, and if a lower function value has been  
c        obtained but the decrease is not sufficient.
c
         if (stage1 .and. f .le. fx .and. f .gt. ftest1) then
c
c           Define the modified function and derivative values.
c
            fm = f - stp*dgtest
            fxm = fx - stx*dgtest
            fym = fy - sty*dgtest
            dgm = dg - dgtest
            dgxm = dgx - dgtest
            dgym = dgy - dgtest
c 
c           Call cstep to update the interval of uncertainty 
c           and to compute the new step.
c
            call cstep(stx,fxm,dgxm,sty,fym,dgym,stp,fm,dgm,
     *                 brackt,stmin,stmax,infoc)
c
c           Reset the function and gradient values for f.
c
            fx = fxm + stx*dgtest
            fy = fym + sty*dgtest
            dgx = dgxm + dgtest
            dgy = dgym + dgtest
         else
c 
c           Call cstep to update the interval of uncertainty 
c           and to compute the new step.
c
            call cstep(stx,fx,dgx,sty,fy,dgy,stp,f,dg,
     *                 brackt,stmin,stmax,infoc)
            end if
c
c        Force a sufficient decrease in the size of the
c        interval of uncertainty.
c
         if (brackt) then
            if (abs(sty-stx) .ge. p66*width1) 
     *         stp = stx + p5*(sty - stx)
            width1 = width
            width = abs(sty-stx)
            end if
c
c        End of iteration.
c
         go to 30
c
c     Last card of subroutine cvsrch.
c
      end
      subroutine cstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
     *                 stpmin,stpmax,info)
      integer info
      double precision stx,fx,dx,sty,fy,dy,stp,fp,dp,stpmin,stpmax
      logical brackt,bound
c     **********
c
c     Subroutine cstep
c
c     The purpose of cstep is to compute a safeguarded step for
c     a linesearch and to update an interval of uncertainty for
c     a minimizer of the function.
c
c     The parameter stx contains the step with the least function
c     value. The parameter stp contains the current step. It is
c     assumed that the derivative at stx is negative in the
c     direction of the step. If brackt is set true then a
c     minimizer has been bracketed in an interval of uncertainty
c     with endpoints stx and sty.
c
c     The subroutine statement is
c
c       subroutine cstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
c                        stpmin,stpmax,info)
c 
c     where
c
c       stx, fx, and dx are variables which specify the step,
c         the function, and the derivative at the best step obtained
c         so far. The derivative must be negative in the direction
c         of the step, that is, dx and stp-stx must have opposite 
c         signs. On output these parameters are updated appropriately.
c
c       sty, fy, and dy are variables which specify the step,
c         the function, and the derivative at the other endpoint of
c         the interval of uncertainty. On output these parameters are 
c         updated appropriately.
c
c       stp, fp, and dp are variables which specify the step,
c         the function, and the derivative at the current step.
c         If brackt is set true then on input stp must be
c         between stx and sty. On output stp is set to the new step.
c
c       brackt is a logical variable which specifies if a minimizer
c         has been bracketed. If the minimizer has not been bracketed
c         then on input brackt must be set false. If the minimizer
c         is bracketed then on output brackt is set true.
c
c       stpmin and stpmax are input variables which specify lower 
c         and upper bounds for the step.
c
c       info is an integer output variable set as follows:
c         If info = 1,2,3,4,5, then the step has been computed
c         according to one of the five cases below. Otherwise
c         info = 0, and this indicates improper input parameters.
c
c     Subprograms called
c
c       FORTRAN-supplied ... abs,max,min,sqrt
c                        ... dble
c
c     Argonne National Laboratory. MINPACK Project. June 1983
c     Jorge J. More', David J. Thuente
c
c     **********
      double precision gamma,p,p66,q,r,s,sgnd,stpc,stpf,stpq,theta
      data p66 /0.66/
      info = 0
c
c     Check the input parameters for errors.
c
      if ((brackt .and. (stp .le. min(stx,sty) .or.
     *    stp .ge. max(stx,sty))) .or.
     *    dx*(stp-stx) .ge. 0.0 .or. stpmax .lt. stpmin) return
c
c     Determine if the derivatives have opposite sign.
c
      sgnd = dp*(dx/abs(dx))
c
c     First case. A higher function value.
c     The minimum is bracketed. If the cubic step is closer
c     to stx than the quadratic step, the cubic step is taken,
c     else the average of the cubic and quadratic steps is taken.
c
      if (fp .gt. fx) then
         info = 1
         bound = .true.
         theta = 3*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
         gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
         if (stp .lt. stx) gamma = -gamma
         p = (gamma - dx) + theta
         q = ((gamma - dx) + gamma) + dp
         r = p/q
         stpc = stx + r*(stp - stx)
         stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/2)*(stp - stx)
         if (abs(stpc-stx) .lt. abs(stpq-stx)) then
            stpf = stpc
         else
           stpf = stpc + (stpq - stpc)/2
           end if
         brackt = .true.
c
c     Second case. A lower function value and derivatives of
c     opposite sign. The minimum is bracketed. If the cubic
c     step is closer to stx than the quadratic (secant) step, 
c     the cubic step is taken, else the quadratic step is taken.
c
      else if (sgnd .lt. 0.0) then
         info = 2
         bound = .false.
         theta = 3*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
         gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
         if (stp .gt. stx) gamma = -gamma
         p = (gamma - dp) + theta
         q = ((gamma - dp) + gamma) + dx
         r = p/q
         stpc = stp + r*(stx - stp)
         stpq = stp + (dp/(dp-dx))*(stx - stp)
         if (abs(stpc-stp) .gt. abs(stpq-stp)) then
            stpf = stpc
         else
            stpf = stpq
            end if
         brackt = .true.
c
c     Third case. A lower function value, derivatives of the
c     same sign, and the magnitude of the derivative decreases.
c     The cubic step is only used if the cubic tends to infinity 
c     in the direction of the step or if the minimum of the cubic
c     is beyond stp. Otherwise the cubic step is defined to be 
c     either stpmin or stpmax. The quadratic (secant) step is also 
c     computed and if the minimum is bracketed then the the step 
c     closest to stx is taken, else the step farthest away is taken.
c
      else if (abs(dp) .lt. abs(dx)) then
         info = 3
         bound = .true.
         theta = 3*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
c
c        The case gamma = 0 only arises if the cubic does not tend
c        to infinity in the direction of the step.
c
         gamma = s*sqrt(max(0.,(theta/s)**2 - (dx/s)*(dp/s)))
         if (stp .gt. stx) gamma = -gamma
         p = (gamma - dp) + theta
         q = (gamma + (dx - dp)) + gamma
         r = p/q
         if (r .lt. 0.0 .and. gamma .ne. 0.0) then
            stpc = stp + r*(stx - stp)
         else if (stp .gt. stx) then
            stpc = stpmax
         else
            stpc = stpmin
            end if
         stpq = stp + (dp/(dp-dx))*(stx - stp)
         if (brackt) then
            if (abs(stp-stpc) .lt. abs(stp-stpq)) then
               stpf = stpc
            else
               stpf = stpq
               end if
         else
            if (abs(stp-stpc) .gt. abs(stp-stpq)) then
               stpf = stpc
            else
               stpf = stpq
               end if
            end if
c
c     Fourth case. A lower function value, derivatives of the
c     same sign, and the magnitude of the derivative does
c     not decrease. If the minimum is not bracketed, the step
c     is either stpmin or stpmax, else the cubic step is taken.
c
      else
         info = 4
         bound = .false.
         if (brackt) then
            theta = 3*(fp - fy)/(sty - stp) + dy + dp
            s = max(abs(theta),abs(dy),abs(dp))
            gamma = s*sqrt((theta/s)**2 - (dy/s)*(dp/s))
            if (stp .gt. sty) gamma = -gamma
            p = (gamma - dp) + theta
            q = ((gamma - dp) + gamma) + dy
            r = p/q
            stpc = stp + r*(sty - stp)
            stpf = stpc
         else if (stp .gt. stx) then
            stpf = stpmax
         else
            stpf = stpmin
            end if
         end if
c
c     Update the interval of uncertainty. This update does not
c     depend on the new step or the case analysis above.
c
      if (fp .gt. fx) then
         sty = stp
         fy = fp
         dy = dp
      else
         if (sgnd .lt. 0.0) then
            sty = stx
            fy = fx
            dy = dx
            end if
         stx = stp
         fx = fp
         dx = dp
         end if
c
c     Compute the new step and safeguard it.
c
      stpf = min(stpmax,stpf)
      stpf = max(stpmin,stpf)
      stp = stpf
      if (brackt .and. bound) then
         if (sty .gt. stx) then
            stp = min(stx+p66*(sty-stx),stp)
         else
            stp = max(stx+p66*(sty-stx),stp)
            end if
         end if
      return
c
c     Last card of subroutine cstep.
c
      end


