      module switches
         integer, parameter :: ipsr     = 2      ! real: 1=decoherence, 2=broadening
         integer, parameter :: ipsv     = 1      ! virt: 1=decoherence, 2=broadening
         integer, parameter :: qqq      = 1      ! scenario: 1=sudakov, 2=no quench, 3=normal
         logical, parameter :: dla      = .false. ! checking DLA formulae
         integer, parameter :: lin      = 2      ! diff step: 1=newton step, 2=rk4
         integer, parameter :: integral = 1      ! integral: 1=trapezoidal, 2=Simpsons
      end module switches

      module scales
         real(8), parameter :: q0=0.25d0
         real(8), parameter :: qhat0=1.5d0 ! Gev^2/fm
         real(8), parameter :: length0=2.d0 ! fm
         real(8), parameter :: gevfm=5.d0
         real(8), parameter :: emin=1.d0
         real(8), parameter :: pi=3.141593d0
         real(8), parameter :: nspec=6.d0
         real(8), parameter :: ca=3.d0, cf=1.333333d0
         ! real(8), parameter :: beta=4.5d0 ! (11/3 Nc - 4/3 nf TR)/2
         real(8), parameter :: b0=0.7161972d0 ! (33-2nf)/(12pi)
         real(8), parameter :: lqcd=0.15d0
         real(8), parameter :: abarmed=0.3d0
         real(8) :: qhat,length,thetac,omr,omc
         real(8) :: ejet,rjet
         real(8) :: nu_int
         save
         contains

            subroutine setscales(ee,rr)
            real(8) :: ee,rr
            ejet = ee
            rjet = rr
            qhat = qhat0/gevfm
            length = length0*gevfm
            thetac = (12.d0/qhat/length**3.d0)**0.5d0
            omc = qhat*length**2.d0
            omr = dsqrt(qhat*length)/rjet
            ! print*,"thetac = ",thetac," omd = ",qhat*length**2.d0
            end subroutine setscales

            subroutine quenching(qfac)
            real(8), intent(out) :: qfac
            real(8) :: qfac1
            qfac = qfac1(ejet)
            end subroutine quenching

            double precision function alphas(xx)
            use switches
            real(8) :: xx
            if(dla) then
               alphas = 0.3d0
            else 
               alphas = 0.5d0/b0
     +            / (dlog(xx/q0) + 0.d0*dlog(q0/lqcd))
            endif
            return
            end      

            double precision function pgg(xx)
            use switches
            real(8) :: xx
            if(dla) then
               pgg = ca/xx
            else
               pgg = ca*(1.d0-xx*(1.d0-xx))**2.d0/xx/(1.d0-xx)
            endif
            return
            end
      
            double precision function pgq(xx)
            use switches
            real(8) :: xx
            if(dla) then
               pgq = cf/xx
            else
               pgq = cf*(1.d0+(1.d0-xx)**2.d0)/xx
            endif
            return
            end
            
       end module scales
       
       module grids
         integer, parameter :: n=2000,m=50
         real(8), dimension(n+1) :: pval,deltap
         real(8), dimension(m+1) :: yval,deltay
         real(8), dimension(n+1,m+1) :: dist
         real(8), dimension(n+1,m+1,n+1) :: ww
         save
         contains
         
            subroutine readic
            implicit none
            integer :: i,k
                        
            do i = 1,n+1
               do k = 1,m+1
                  dist(i,k) = 0.d0
               enddo
            enddo
            do i = 1,n+1
               dist(i,1) = 1.d0
            enddo
            end subroutine readic

            subroutine init
            ! setting up grids for the distribution
            use scales
            real(8) :: zz,kt,theta
            real(8) :: pmin,pmax,parg,test,q1

            ! q1 = 1.01d0*q0
            pmin = emin
            pmax = ejet
            ! linear grid      
            ! do i = 1,n+1
            !    deltap(i) = (pmax - pmin) / real(n)
            !    pval(i) = pmin + real(i-1)*deltap(i)
            ! enddo
            ! symmetric logarithmic grid
   !          do i = 1,n+1
   !             parg = real(2*(i-1)-n) / real(n) * dlog(pmax/pmin)
   !             pval(i) = dexp(parg) / (1.d0 + dexp(parg))
   !   &         * (pmax + pmin)
   !             deltap(i) = 2.d0/real(n) 
   !   &         * dlog(pmax/pmin) * pval(i)
   !   &         / (1.d0 + dexp(parg))
   !          enddo
            ! asymmetric logarithmic grid
            do i = 1,n+1
               pval(i) = pmin*(pmax/pmin)**(real(i-1)/real(n))
               deltap(i) = pval(i)
     +            *((pmax/pmin)**(1.d0/real(n)) - 1.d0)
            enddo
            ! linear grid for y :: [ylow,yhigh]
            ymin = dlog(rjet/thetac) ! initial condition (minimal angle)
            ymax =  0.d0             ! jet radius (maximal angle)
            do k = 1,m+1
               deltay(k) = -1.d0*(ymax - ymin) / real(m)
               yval(k) = ymin - real(k-1)*deltay(k)
            enddo
            
            call readic
            do i=1,n+1
               do k=1,m+1
                  theta = rjet*dexp(-yval(k))
                  do j=1,n+1
                     zz = pval(j)/pval(i)
                     if(zz.le.1.d0.and.zz.ge.0.d0) then
                        kt = zz*(1.d0-zz)*pval(i)*theta
                        if(kt.gt.1.05d0*q0) then
                           ww(i,k,j) = alphas(kt)*pgg(zz)/pval(i)/pi
                        else
                           ww(i,k,j) = 0.d0
                        endif
                     else
                        ww(i,k,j) = 0.d0
                     endif
                  enddo
               enddo
            enddo
            end subroutine init
      end module grids

      program collimator_gluon_eloss
c     ****************************************************
c     Program to solve the non-linear evolution equation
c     for the collimator function.
c     Computing for set of Rjet and pTjet.
c     ****************************************************
      use grids
      use scales
      use switches
      implicit none
      real(8) :: pjet,pmax,pmin
      real(8) :: thetajet,thetamax,thetamin,theta
      real(8) :: qf
      real(8), dimension(100) :: pp,raa
      integer :: i,j,k
      integer :: nplot1,nplot2
         
      ! Plotting parameters
      pmax=10000.d0
      pmin=10.d0
      thetamin=0.3d0
      thetamax=0.8d0
      nplot1=6
      nplot2=30

      print*,"Quenching factor"
      do i = 1,nplot2
         pjet = pmin*(pmax/pmin)**(real(i-1)/real(nplot2-1))
         pp(i) = pjet
         call setscales(pp(i),0.4d0)
         call quenching(qf)
         print*,pjet,qf
      enddo
      do j = 1,nplot1
         thetajet = thetamin+(thetamax-thetamin)
     &      *real(j-1)/real(nplot1-1)
         print*,"RAA at R=",thetajet
         do i = 1,nplot2
            call setscales(pp(i),thetajet)
            ! call setscales(10.d0,thetajet)
            call init
c     **********************************************************************
ccc   Evolution in y
            do k = 2,m+1
               call rk4_step(k)
            enddo
c     **********************************************************************
         print*,pp(i),dist(n+1,m+1)
         ! print*,thetajet,",",dist(n+1,m+1)
         ! do k = 1,n+1,20
         !    print*,dist(k,m+1)
         ! enddo   
         enddo
         ! do i = 1,nplot2
         !    print*,pp(i),raa(i)
         ! enddo
      enddo
      end
      

      double precision function qfac1(pt)
      use scales
      use switches
      real(8) :: dgauss1,pt
      external emult
      if(qqq.eq.1) then
         qfac1 = 0.d0
      elseif(qqq.eq.2) then
         qfac1 = 1.d0
      elseif(qqq.eq.3) then
         nu_int = nspec/pt
         res = dgauss1(emult,0.d0,omr,1.d-5)
         qfac1 = dexp(-nu_int*res)
      else
         print*,"Uknkown quench option!"
         stop
      endif
      return
      end

      double precision function emult(xx)
      use scales
      real(8) :: xx
      double precision mult
      emult = dexp(-nu_int*xx)*mult(xx)
      return
      end

      double precision function mult(xx)
      use scales
      real(8) :: xx
      mult = 2.d0*abarmed*(dsqrt(omc/xx)
     &   + dlog(2.d0)*dlog(2.d0*xx/omc)
     &   - 1.44136d0)
      return
      end

      subroutine rk4_step(k)
ccc   4th order Runge-Kutta
      use grids
      use switches
      implicit none
      real(8), dimension(n+1) :: y,y1
      real(8), dimension(n+1) :: dydx,dydx1,dydx2,dydx3
      real(8) :: dy,dy2,dy6
      real(8) :: aa1
      integer :: i,k
                       
      dy  = deltay(k)
      dy2 = deltay(k) / 2.d0
      dy6 = deltay(k) / 6.d0

      do i = 1,n+1
            y(i) = dist(i,k-1)
      enddo

      ! print*,"going into deriv1:",k
      call deriv2(k,dydx,y)      ! k1
      if(lin.eq.1) then
         do i = 1,n+1
           dist(i,k) = y(i) + dy*dydx(i)
         enddo
      elseif(lin.eq.2) then
         do i = 1,n+1
            y1(i) = y(i) + dy2 * dydx(i)
         enddo
      ! print*,"going into deriv2:",k
         call deriv2(k,dydx1,y1)    ! k2
         do i = 1,n+1
            y1(i) = y(i) + dy2 * dydx1(i)
         enddo
      ! print*,"going into deriv3:",k
         call deriv2(k,dydx2,y1)    ! k3
         do i = 1,n+1
            y1(i) = y(i) + dy * dydx2(i)
            dydx3(i) = dydx2(i) + dydx1(i)
         enddo
      ! print*,"going into deriv4:",k
         call deriv2(k,dydx1,y1)    ! k4
         do i = 1,n+1
            dist(i,k) = y(i) + dy6
     &         * (dydx(i) + dydx1(i) + 2.d0 * dydx3(i))
         enddo
      else
         print*,"Unkown option for lin!"
      endif
      end
      
      subroutine deriv2(k1,dydx1,y1)
      ! implementing composite Simpson's rule for irregulary spaced data
      ! https://en.wikipedia.org/wiki/Simpson%27s_rule
      ! (note that there are typos in the alpha and beta coefficients,
      ! should be order h^3 in all terms in the numerator!)
      ! in case of an ODD number of subintervals: trapezoidal for the last interval
      use grids
      use scales
      use switches
      implicit none
      real(8), dimension(n+1), intent(in) :: y1
      real(8), dimension(n+1), intent(out) :: dydx1
      real(8), dimension(n+1) :: dh,dl
      real(8) :: al,be,nu,dda,ddb,ddn
      real(8) :: mish,misl
      real(8) :: qfac
      integer, dimension(n+1) :: ilow,ihigh
      integer :: i,j,k,k1,imax,imin
      integer :: ilim,deli,jj,ipsr1,ipsv1

      call quenching(qfac)
      qfac = qfac**2.d0

      do i=1,n+1
         dydx1(i) = 0.d0
      enddo
      k=k1

ccc   Real term
      ! Finding the limits for integral contribution
      ipsr1=ipsr
      call indices(k,ipsr1,ilim,ilow,ihigh,dh,dl)

      do i=ilim,n+1
         imax=ihigh(i)
         imin=ilow(i)
         deli=imax-imin ! number of subintervals
         if(deli.le.0) then
            continue
         elseif(deli.eq.1) then ! 1 interval (trapezoidal)
            dda = ww(i,k,imin)*y1(imin)*y1(i-imin)*qfac
            ddb = ww(i,k,imax)*y1(imax)*y1(i-imax)*qfac
            dydx1(i) = dydx1(i) 
     &      + 0.5d0*(pval(imax)-pval(imin))*(dda+ddb)
         elseif(deli.ge.2) then ! 2 or more subintervals
            if(integral.eq.1) then ! trapezoidal
               do j=1,deli
                  jj = j+imin-1
                  dda = ww(i,k,jj)*y1(jj)*y1(i-jj)*qfac
                  ddb = ww(i,k,jj+1)*y1(jj+1)*y1(i-jj-1)*qfac
                  if(j.eq.1) then
                     dydx1(i) = dydx1(i) 
     &                  + 0.5d0*deltap(jj)*(dda+ddb)
                  else
                    dydx1(i) = dydx1(i) 
     &                  + 0.5d0*deltap(jj)*(dda+ddb)
                  endif
               enddo
            elseif(integral.eq.2) then
               if(mod(deli,2).ne.0) then
                  dda = ww(i,k,imax-1)*y1(imax-1)*y1(i-imax+1)*qfac
                  ddb = ww(i,k,imax)*y1(imax)*y1(i-imax)*qfac
                  dydx1(i) = dydx1(i) 
     &               + 0.5d0*(pval(imax)-pval(imax-1))*(dda+ddb)
                  deli = deli - 1
               endif
               ! now we have EVEN number of subinteravls   
               do j=0,deli/2-1
                  jj = 2*j+imin
                  al = 0.166667d0*(deltap(jj+1)+deltap(jj))
     &               * (2.d0*deltap(jj+1)-deltap(jj))/deltap(jj+1)
                  be = 0.166667d0*(deltap(jj+1)+deltap(jj))**3.d0
     &               / deltap(jj+1) / deltap(jj)
                  nu = 0.166667d0*(deltap(jj+1)+deltap(jj))
     &               * (2.d0*deltap(jj)-deltap(jj+1))/deltap(jj)
                  ddn = ww(i,k,jj)*y1(jj)*y1(i-jj)*qfac
                  ddb = ww(i,k,jj+1)*y1(jj+1)*y1(i-jj-1)*qfac
                  dda = ww(i,j,jj+2)*y1(jj+2)*y1(i-jj-2)*qfac
                  dydx1(i) = dydx1(i) + al*dda + be*ddb + nu*ddn
               enddo
            else
               print*,"Unkown integration option!"
               stop
            endif   
         endif
         mish = dh(i)*ww(i,k,imax)*y1(imax)*y1(i-imax)*qfac
         misl = dl(i)*ww(i,k,imin)*y1(imin)*y1(i-imin)*qfac
         dydx1(i) = dydx1(i) + mish + misl
      enddo
 
ccc   Virtual term
      ! Finding the limits for integral contribution
      ipsv1=ipsv
      if(ipsv.ne.ipsr) then
         call indices(k,ipsv1,ilim,ilow,ihigh,dh,dl)
      endif

      do i=ilim,n+1
         imax=ihigh(i)
         imin=ilow(i)
         deli=imax-imin ! number of subintervals
         if(deli.le.0) then
            continue
         elseif(deli.eq.1) then ! 1 interval (trapezoidal)
            dda = -ww(i,k,imin)*y1(i)
            ddb = -ww(i,k,imax)*y1(i)
            dydx1(i) = dydx1(i) 
     &      + 0.5d0*(pval(imax)-pval(imin))*(dda+ddb)      
         elseif(deli.ge.2) then ! 2 or more subintervals
            if(integral.eq.1) then ! trapezoidal
               do j=1,deli
                  jj = j+imin-1
                  dda = -ww(i,k,jj)*y1(i)
                  ddb = -ww(i,k,jj+1)*y1(i)
                  dydx1(i) = dydx1(i) + 0.5d0*deltap(jj)*(dda+ddb)      
               enddo
            elseif(integral.eq.2) then ! Simpson
               if(mod(deli,2).ne.0) then
                  dda = -ww(i,k,imax-1)*y1(i)
                  ddb = -ww(i,k,imax)*y1(i)
                  dydx1(i) = dydx1(i)
     &               + 0.5d0*(pval(imax)-pval(imax-1))*(dda+ddb)
                  deli = deli - 1
               endif
c           now we have EVEN number of subinteravls   
               do j=0,deli/2-1
                  jj = 2*j+imin
                  al = 0.166667d0*(deltap(jj+1)+deltap(jj))
     &               * (2.d0*deltap(jj+1)-deltap(jj))/deltap(jj+1)
                  be = 0.166667d0*(deltap(jj+1)+deltap(jj))**3.d0
     &               / deltap(jj+1) / deltap(jj)
                  nu = 0.166667d0*(deltap(jj+1)+deltap(jj))
     &               * (2.d0*deltap(jj)-deltap(jj+1))/deltap(jj)
                  ddn = -ww(i,k,jj)*y1(i)
                  ddb = -ww(i,k,jj+1)*y1(i)
                  dda = -ww(i,j,jj+2)*y1(i)
                  dydx1(i) = dydx1(i) + al*dda + be*ddb + nu*ddn
               enddo
            else
               print*,"Unkown integration option!"
               stop
            endif   
         endif
         mish = -dh(i)*y1(i)*ww(i,k,imax)
         misl = -dl(i)*y1(i)*ww(i,k,imin)
         dydx1(i) = dydx1(i) + mish + misl
      enddo
      end

      subroutine indices(k,ips,ilim,ilow,ihigh,dh,dl)
      use grids
      use scales
      implicit none
      real(8), dimension(n+1), intent(out) :: dh,dl
      real(8) :: aa1,aa2,om,theta,omplu,ommin
      integer, dimension(n+1), intent(out) :: ilow,ihigh
      integer :: i,j,k,ips,ilim

      theta = rjet*dexp(-yval(k))
      ! aa2 = 2.d0/length/theta**2.d0 ! length
      if(ips.eq.1) then ! decoherence
         aa1 = (0.66667d0*qhat/theta**4.d0)**0.33333d0
      elseif(ips.eq.2) then ! broadening
         ! aa1 = (qhat*length)**0.5d0/theta
         if(theta.lt.(thetac*rjet**3.d0)**0.25d0) then
            aa1 = (0.66667d0*qhat/theta**4.d0)**0.33333d0
         else
            aa1 = (0.333d0*qhat*length)**0.5d0/rjet ! factor 1/3
            ! aa1 = (qhat*length)**0.5d0/rjet
         endif
      else
         print*,"Unknown ips..."
         stop
      endif

      i=1
      do while(pval(i).lt.4.d0*aa1
     &   .and. i.le.n+1)
         i=i+1
         ilow(i)=0
         ihigh(i)=0
         dh(i)=0.d0
         dl(i)=0.d0
      enddo
      ilim = i

      do i=ilim,n+1
ccc   finding the upper and lower limit on the integral
         om = pval(i)
         omplu = 0.d0
         ommin = 0.d0
         if(om > 4.d0*aa1) then
            omplu = om*0.5d0*(1.d0+dsqrt(1.d0-4.d0*aa1/om))
            ommin = om*0.5d0*(1.d0-dsqrt(1.d0-4.d0*aa1/om))
         else 
            stop
         endif
         j=1 ! starting at omega'=0
         do while(pval(j).lt.ommin.and.j.le.n+1)
            j=j+1
         enddo
         ilow(i)=j
c         print*,"  ",imin,"  ",pval(imin),"   ",ommin
         
         j=i ! starting right below omega
         do while(pval(j).gt.omplu.and.j.ge.1)
            j=j-1
         enddo
         ihigh(i)=j

         dh(i) = omplu-pval(ihigh(i))
         dl(i) = pval(ilow(i))-ommin
         ! print*,"min ",ommin,pval(ilow(i))
         ! print*,"plu ",omplu,pval(ihigh(i))
      enddo
      end


      DOUBLE PRECISION FUNCTION DGAUSS1(F,A,B,EPS)
      DOUBLE PRECISION F,A,B,EPS
      DOUBLE PRECISION W(12),X(12),AA,BB,C1,C2,U,S8,S16,CONST
      LOGICAL MFLAG,RFLAG
      EXTERNAL F
C
C     ******************************************************************
C
C     ADAPTIVE DOUBLE PRECISION GAUSSIAN QUADRATURE.
C
C     DGAUSS1 IS SET EQUAL TO THE APPROXIMATE VALUE OF THE INTEGRAL OF
C     THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY PARAMETER
C     EPS.
C
C     ******************************************************************
C
      DATA W / 0.10122 85362 90376 259D0,
     1         0.22238 10344 53374 471D0,
     2         0.31370 66458 77887 287D0,
     3         0.36268 37833 78361 983D0,
     4         0.27152 45941 17540 949D-1,
     5         0.62253 52393 86478 929D-1,
     6         0.95158 51168 24927 848D-1,
     7         0.12462 89712 55533 872D0,
     8         0.14959 59888 16576 732D0,
     9         0.16915 65193 95002 538D0,
     A         0.18260 34150 44923 589D0,
     B         0.18945 06104 55068 496D0/

      DATA X / 0.96028 98564 97536 232D0,
     1         0.79666 64774 13626 740D0,
     2         0.52553 24099 16328 986D0,
     3         0.18343 46424 95649 805D0,
     4         0.98940 09349 91649 933D0,
     5         0.94457 50230 73232 576D0,
     6         0.86563 12023 87831 744D0,
     7         0.75540 44083 55003 034D0,
     8         0.61787 62444 02643 748D0,
     9         0.45801 67776 57227 386D0,
     A         0.28160 35507 79258 913D0,
     B         0.95012 50983 76374 402D-1/
C
C     ******************************************************************
C
C  START.
      DGAUSS1=0.0D0
      IF(B.EQ.A) RETURN
      CONST=0.005D0/(B-A)
      BB=A
C
C  COMPUTATIONAL LOOP.
    1 AA=BB
      BB=B
    2    C1=0.5D0*(BB+AA)
         C2=0.5D0*(BB-AA)
         S8=0.0D0
         DO 3 I=1,4
            U=C2*X(I)
            S8=S8+W(I)*(F(C1+U)+F(C1-U))
    3    CONTINUE
         S8=C2*S8
         S16=0.0D0
         DO 4 I=5,12
            U=C2*X(I)
            S16=S16+W(I)*(F(C1+U)+F(C1-U))
    4    CONTINUE
         S16=C2*S16
         IF( ABS(S16-S8) .LE. EPS*(1.+ABS(S16)) ) GO TO 5
         BB=C1
         IF( 1.D0+ABS(CONST*C2) .NE. 1.D0) GO TO 2
      DGAUSS1=0.0D0
      ! CALL KERMTR('D103.1',LGFILE,MFLAG,RFLAG)
      ! IF(MFLAG) THEN
      !    IF(LGFILE.EQ.0) THEN
      ! WRITE(*,6)
      !    ELSE
      !       WRITE(LGFILE,6)
      !    ENDIF
      ! ENDIF
      ! IF(.NOT. RFLAG) CALL ABEND
      RETURN
    5 DGAUSS1=DGAUSS1+S16
      IF(BB.NE.B) GO TO 1
      RETURN
C
    6 FORMAT( 4X, 'FUNCTION DGAUSS1 ... TOO HIGH ACCURACY REQUIRED')
      END


! ccc   old code: assuming equidistant intervals
!             dd = dist(imin)*dist(i-imin)*qfac - dist(i)
!             xdist = deltap(imin)*ww(i,k,imin)*dd /3.d0               
!             do j=imin+1,imax-1
!                   dd = dist(j)*dist(i-j)*qfac - dist(i)
!                   xdist = xdist
!      &               + 2.d0 * (1.d0 + mod(j+1,2))/3.d0
!      &               * deltap(j)*ww(i,k,j)*dd
!                enddo
!                dd = dist(imax)*dist(i-imax)*qfac - dist(i)
!                xdist = xdist + deltap(imax)*ww(i,k,imax)*dd /3.d0
!             elseif(mod(deli,2).eq.0) then ! even # of intermediate
! c              initially do a trapezoidal step               
!                dd = dist(imin)*dist(i-imin)*qfac - dist(i)
!                xdist = deltap(imin)*ww(i,k,imin)*dd /2.d0
!                dd = dist(imin+1)*dist(i-imin-1)*qfac - dist(i)
!                xdist = xdist + deltap(imin+1)*ww(i,k,imin+1)*dd /2.d0
! c              continue with Simpson      
!                xdist = xdist + deltap(imin+1)*ww(i,k,imin+1)*dd /3.d0
!                do j=imin+2,imax-1
!                   dd = dist(j)*dist(i-j)*qfac - dist(i)
!                   xdist = xdist
!      &               + 2.d0 * (1.d0 + mod(j+1,2))/3.d0
!      &               * deltap(j)*ww(i,k,j)*dd
!                enddo
!                dd = dist(imax)*dist(i-imax)*qfac - dist(i)
!                xdist = xdist + deltap(imax)*ww(i,k,imax)*dd /3.d0
! ccc


!       subroutine deriv(k,dydx,y1)
!       implicit none
!       double precision dydx,dist,ww
!       double precision al,be,nu,dda,ddb,ddn
!       double precision om,thet,aa1
!       double precision omplu,ommin
!       double precision pval,yval,deltap,deltay
!       double precision qhat,length,rjet,ejet
!       double precision qfac
!       integer i,j,n,m,k
!       integer ilim,imin,imax,deli,dint,jj
!       integer integral

!       parameter(n=2000,m=200)
            
!       dimension dydx(1:n+1),y1(1:n+1)
!       dimension pval(1:n+1),yval(1:m+1)
!       dimension deltap(1:n+1),deltay(1:m+1)
!       dimension ww(1:n+1,1:m+1,1:n+1)
      
!       common/steps/ deltap,deltay
!       common/grid/ pval,yval
!       common/weights/ ww
!       common/scales/ qhat,length,rjet,ejet
!       common/quench/ qfac

!       data integral/2/

! ccc   Finding the limit for integral contribution
!       thet = rjet*dexp(-yval(k))
!       aa1 = (0.66667d0*qhat/thet**4.d0)**0.33333d0
! c      print*,"At angle: ",thet," aa1 = ",aa1
!       i=1
!       do while(pval(i).lt.4.d0*aa1.and.i.lt.n+1)
!          dydx(i)=0.d0
!          i=i+1
!       enddo
!       ilim=i
! c      print*,ilim," ",pval(ilim)," ",4.d0*aa1

!       do i=ilim,n+1
! ccc   finding the upper and lower limit on the integral
!          om = pval(i)
! c         print*,i," ",om
!          omplu = 0.d0
!          ommin = 0.d0
!          if(om > 4.d0*aa1) then
!             omplu = om*0.5d0*(1.d0+dsqrt(1.d0-4.d0*aa1/om))
!             ommin = om*0.5d0*(1.d0-dsqrt(1.d0-4.d0*aa1/om))
!          else 
!             stop
!          endif
!          j=1 ! starting at omega'=0
!          do while (pval(j).lt.ommin.and.j.lt.n+1)
!             j=j+1
!          enddo
!          imin=j
! c         print*,"  ",imin,"  ",pval(imin),"   ",ommin

!          j=i ! starting right below omega
!          do while (pval(j).gt.omplu.and.j.gt.0)
!             j=j-1
!          enddo
!          imax=j
! c         print*,"  ",imax,"  ",pval(imax),"   ",omplu

!          deli=imax-imin ! number of subintervals
!          ! print*,"   numer of subintervals:",deli

!          if(deli.le.0) then
!             ! print*,"Error in computing the integral!"
!             ! print*,"  ",imin,"  ",pval(imin),"   ",ommin
!             ! print*,"  ",imax,"  ",pval(imax),"   ",omplu
!             ! stop
!          ! elseif(deli.eq.0) then ! no interval
!             dydx(i) = 0.d0
!          elseif(deli.eq.1) then ! 1 interval (trapezoidal)
!             dda = ww(i,k,imin)
!      &          * (dist(imin)*dist(i-imin)*qfac - dist(i))
!             ddb = ww(i,k,imax)
!      &          * (dist(imax)*dist(i-imax)*qfac - dist(i))
!             dydx(i) = 0.5d0*(pval(imax)-pval(imin))
!      &          * (dda + ddb)      
! c            print*,"trap:",dda,ddb
!          elseif(deli.ge.2) then ! 2 or more subintervals
!             dydx(i) = 0.d0
!             if(integral.eq.1) then ! trapezoidal
!                do j=1,deli
!                   jj = j+imin-1
!                   dda = ww(i,k,jj)
!      &               * (dist(jj)*dist(i-jj)*qfac - dist(i))
!                   ddb = ww(i,k,jj+1)
!      &               * (dist(jj+1)*dist(i-jj-1)*qfac - dist(i))
!                   dydx(i) = dydx(i) + 0.5d0*deltap(jj) * (dda + ddb)      
!                enddo
!             elseif(integral.eq.2) then
! c           implementing composite Simpson's rule for irregulary spaced data
! c           https://en.wikipedia.org/wiki/Simpson%27s_rule
! c           (note that there are typos in the alpha and beta coefficients,
! c           should be order h^3 in all terms in the numerator!)
! c           in case of an ODD number of subintervals: treat last interval in a special way
!                if(mod(deli,2).ne.0) then
!    !             al = (2.d0*deltap(imax-1)**2.d0
!    !   &            + 3.d0*deltap(imax-1)*deltap(imax-2))
!    !   &            / (6.d0*(deltap(imax-2)+deltap(imax-1)))
!    !             be = (deltap(imax-1)**2.d0
!    !   &            + 3.d0*deltap(imax-1)*deltap(imax-2))
!    !   &            / (6.d0*deltap(imax-2))
!    !             nu = - deltap(imax-1)**3.d0
!    !   &            / (6.d0*deltap(imax-2)
!    !   &            * (deltap(imax-2)+deltap(imax-1)))
!    !             al = 0.416667d0*deltap(imax)
!    !             be = 0.666667d0*deltap(imax)
!    !             nu = -0.083333d0*deltap(imax)

!    !             dda = ww(i,k,imax)
!    !   &             * (dist(imax)*dist(i-imax)*qfac - dist(i))
!    !             ddb = ww(i,k,imax-1)
!    !   &             * (dist(imax-1)*dist(i-imax+1)*qfac - dist(i))
!    !             ddn = ww(i,k,imax-2)
!    !   &             * (dist(imax-2)*dist(i-imax+2)*qfac - dist(i))
!    !             dydx(i) = al*dda + be*ddb + nu*ddn

!                   dda = ww(i,k,imax-1)
!      &               * (dist(imax-1)*dist(i-imax+1)*qfac - dist(i))
!                   ddb = ww(i,k,imax)
!      &               * (dist(imax)*dist(i-imax)*qfac - dist(i))
!                   dydx(i) = 0.5d0*(pval(imax)-pval(imax-1))
!      &               * (dda + ddb)      
          
!                   deli = deli - 1
!                endif
! c           if(mod(deli,2).ne.0) then ! odd # of intermediate (Simpson)
! c           now we have EVEN number of subinteravls   
!                do j=0,deli/2-1
!                   jj = 2*j+imin
!    !             al = (2.d0*deltap(jj+1)**3.d0 
!    !   &            - deltap(jj)**3.d0
!    !   &            + 3.d0*deltap(jj)*deltap(jj+1))
!    !   &            / (6.d0*deltap(jj+1)
!    !   &            * (deltap(jj+1)+deltap(jj)))
!    !             be = (deltap(jj+1)**3.d0
!    !   &            + deltap(jj)**3.d0
!    !   &            + 3.d0*deltap(jj+1)*deltap(jj)
!    !   &            * (deltap(jj+1)+deltap(jj)))
!    !   &            / (6.d0*deltap(jj+1)*deltap(jj))
!    !             nu = (2.d0*deltap(jj)**3.d0
!    !   &            - deltap(jj+1)**3.d0
!    !   &            + 3.d0*deltap(jj+1)*deltap(jj)**2.d0)
!    !   &            / (6.d0*deltap(jj)
!    !   &            * (deltap(jj+1)+deltap(jj)))
!                   al = 0.166667d0*(deltap(jj+1)+deltap(jj))
!      &               * (2.d0*deltap(jj+1)-deltap(jj))/deltap(jj+1)
!                   be = 0.166667d0*(deltap(jj+1)+deltap(jj))**3.d0
!      &               / deltap(jj+1) / deltap(jj)
!                   nu = 0.166667d0*(deltap(jj+1)+deltap(jj))
!      &               * (2.d0*deltap(jj)-deltap(jj+1))/deltap(jj)
!                ! al = 0.333333d0*deltap(imax)
!                ! be = 1.333333d0*deltap(imax)
!                ! nu = 0.333333d0*deltap(imax)
!                   ddn = ww(i,k,jj)
!      &                * (dist(jj)*dist(i-jj)*qfac - dist(i))
!                   ddb = ww(i,k,jj+1)
!      &                * (dist(jj+1)*dist(i-jj-1)*qfac - dist(i))
!                   dda = ww(i,j,jj+2)
!      &                * (dist(jj+2)*dist(i-jj-2)*qfac - dist(i))
!                   dydx(i) = dydx(i) + al*dda + be*ddb + nu*ddn
!                enddo
!             else
!                print*,"Unkown integration option!"
!                stop
!             endif   
!          endif
!       enddo
!       end

