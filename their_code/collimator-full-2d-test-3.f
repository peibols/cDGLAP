      program    collimator

      implicit none 

      integer, parameter :: n=2000, m=300,n0=20
      real(8), parameter :: Pi=3.14159265D0, Nc=3.D0,nf=5.D0
      real(8) :: CR,CF,CA,beta0,tscale,dtheta,tstep_F,tstep_A,time,taud
      real(8) :: qhat,L,R,omegac,thetac,p0,p1,alpha_F,alpha_A,abar,pmin
      real(8) :: omegas

      real(8), dimension(n)  :: p,jac
      real(8), dimension(n0)  :: y
      real(8), dimension(n,n0)  :: func_F_out,func_A_out,quark,gluon
      real(8), dimension(n0)  :: quark_qw,gluon_qw
      real(8), dimension(n,n) :: split_F,split_A

      
      integer, dimension(n,n) :: iprim
 
      integer, dimension(n0) :: num0

      integer :: grid,i,j,k,diag
      
      character*(18)   arg
      COMMON /input/ p0,p1,qhat,L,alpha_F,thetac,abar

      COMMON /quenching/ quark_qw
      COMMON /radius/ R

      
      CF=(Nc**2-1.D0)/(2*Nc)
      CA=Nc
      beta0=11.D0-2.D0*nf/3.D0

      p0=0.25D0

      p1=5000.D0

      grid=2
      
c     prefactor of the collision integral
      

      tscale=2.D0/beta0
c      tscale=1.D0/CF

c     step in angular integral (time variable)
      
      dtheta=1.D-2
      
c      tstep_F=dtheta*tscale*CF
      tstep_F=dtheta*tscale*CF
c    tstep_A=dtheta*tscale*CA
      tstep_A=dtheta*tscale*CF

      

c     paramters

      qhat=1.5D0/5.D0
      L=15.D0
      R=0.2D0


      alpha_F=CF/Pi*(2*Pi/beta0)
c      alpha_F=0.3
      alpha_A=CA/Pi*(2*Pi/beta0)
c      alpha_A=0.3
      abar=0.3
      
c     relevant scales

c     omega_c

      omegac= qhat * L**2 /12.D0
      omegas=qhat*L**2*(abar)**2.D0
c     theta_c

      thetac= sqrt(1.D0/omegac/L )
      pmin=(2.*qhat/3.)**0.333D0/R**(4.D0/3.D0)



       write(*,*) thetac,pmin,omegas,sqrt(qhat*L)/R,gamma(4.)



***************************************

c      define grid for the integration over the splitting function


       
       call gridvar(p,jac,iprim,split_F,split_A,n,p0,p1,2,n0,y,num0)
        write(*,*) num0(1), split_F(5,5),num0(10)
      write(*,*)   y(5), p(num0(5))
      call qw(gluon_qw,y,n0,5.5D0,alpha_F,qhat,L)
      write(*,*)   y(5), p(num0(5))
      call qw(quark_qw,y,n0,5.5D0,alpha_A,qhat,L)

    
        write(*,*) y(10), gluon_qw(10)

        write(*,*)   y(5), p(num0(5))

********************************************
*   initial condition 
*******************************************
      
        do j=1,n0
           do i=1,n

            quark(i,j)=1.D0
            gluon(i,j)=1.D0

           enddo 
        enddo
********************************************
*   angular evolution  (time)
*******************************************


      time=1.D-10
c      time=thetac



          do k=1,m

         if (mod(k,5) == 1) then
            write(arg,'(F18.10)')  time
         open(10,file='collimator-data-quark/collimator_quark_t_'//arg//
     &         '.dat',status='unknown')
            write(arg,'(F18.10)')  time
         open(20,file='collimator-data-gluon/collimator_gluon_t_'//arg//
     &          '.dat',status='unknown')

              


         
         
         do j=1,n0

            diag=num0(j)
         
             write(10,*)    y(j), quark(diag,j)*quark_qw(j)
c             write(10,*)    y(j), quark(diag,j)
             write(20,*)    y(j), (gluon(diag,j)*gluon_qw(j)
     &      +quark(diag,j)*quark_qw(j))/2.D0
            
c             write(20,*)    y(j), exp(-2*log(R/thetac)
c     &      *(log(y(j)/omegac)+2./3.*log(R/thetac)))

c            write(20,*)    y(j),
c     &   exp(-3./4.*log(R**(4./3.)*y(j)/qhat**0.3333)**2)

          enddo
            
        endif


                  
               
       
                  
          
          call coll(func_F_out,func_A_out,quark,gluon,n,n0,time)

          
      
c        func_F_out=1.D0
c        func_A_out=1.D0

       

         do j=1,n0
           do i=1,n


              quark(i,j) = quark(i,j)+func_F_out(i,j)*tstep_F
              gluon(i,j) = gluon(i,j)+func_A_out(i,j)*tstep_A
            
          enddo
         enddo 

         time=time+dtheta



         close(10)
         close(20)

         
      enddo

********************************************
*   end program
*******************************************

       end



          
***************************************

*  quenching weight   

***************************************

      subroutine  qw(func,x,nstep,power,alpha,qh,length)
         
      integer, parameter :: mstep=10000
      real(8), dimension(nstep), intent(out) :: func 
      real(8), dimension(nstep), intent(in) :: x
      real(8), dimension(nstep) :: integral
      real(8), dimension(nstep,mstep) :: spect,u
      real(8), dimension(mstep) :: w

      real(8), intent(in) :: alpha,length,qh,power
      integer, intent(in) :: nstep
     
      
      real(8) :: Pi,z,y,omegas,omegac,omegaBH,nu,R,temp,a

      COMMON /radius/ R

      Pi=3.14159265D0

      

      omegas=alpha**2*qh*length**2
      omegac=qh*length**2/2
      omegaBH=1.D0
      z=omegaBH/omegas
      z=0.001D0
      y=omegac/omegas

      temp=0.5D0



      
      do j=1,mstep
         
         w(j)=j*5.D0/mstep
         
      enddo

       do i=1,nstep
          do j=1,mstep
             u(i,j)=sqrt(omegac*power/(2*(x(i))*w(j)))
             
             if (u(i,j) > sqrt(omegac/2.D0/omegaBH)) then
                spect(i,j)=0.D0
              else
             spect(i,j) =log ( (cosh(u(i,j)))**2 - (sin(u(i,j)))**2 )
     &    *exp(-(R*x(i)*w(j)/power)**2/qh/length)
             endif
             
         enddo
       enddo
      
      
       
c       spect=log ( (cosh(u)*cos(u))**2 + (sinh(u)*sin(u))**2 )
c      spect= u
      
***   Numerical integration of quenching weights
      
       integral=0.D0

       do i =1,nstep
          do j=1,mstep-1

             integral(i) =  integral(i) + 0.5D0*

             
     &           ( spect(i,j)* (1.D0 - exp(- w(j)))/j
     &           + spect(i,j+1)* (1.D0 - exp(-w(j+1))) / (j+1) )
            

             

            
          enddo
        enddo

c        func=integral

c       func =exp(- alpha * integral/log((2*qh/power*x)**0.25D0/0.25D0))
       func =exp(- 0.3 * integral)
c       a=temp/sqrt(qh*length)*R
c       func=exp(-2.*alpha*sqrt(qh)
c     &   *length*(max(power+0.5*log(x/100.),4.))/x
c     &     *(qh*length)*0.25/R**0.5
c     &     *( gamma(5./4.) - a**0.5+a**2.5/5.))




c      do i=1,nstep
 
  
c        func(i)=exp(-alpha*length*sqrt(Pi*qh*power/x(i)))
 
c      enddo
c      func=1.0D0

      end subroutine

      
***************************************************************       
*     splitting collision term with trapezoidal rule
**************************************************************
      subroutine coll(fun_F_out,fun_A_out,fun_F_in,fun_A_in
     &  ,nstep,nstep0,time0)

      real(8), intent(in) :: time0
      real(8), dimension(nstep,nstep0), intent(in) ::
     &  fun_F_in,fun_A_in
      real(8), dimension(nstep,nstep0), intent(out) ::
     &  fun_F_out,fun_A_out

      real(8), dimension(nstep) :: p,jac,gain
      real(8), dimension(nstep0) :: fun_qw,y0
      real(8), dimension(nstep,nstep) :: split_A,split_F
      real(8), dimension(nstep,nstep) :: collision_A,collision_F
       real(8), dimension(nstep) :: integral_A,integral_F
       integer, dimension(nstep,nstep) :: jprim
       integer, dimension(nstep0) :: num
      real(8), dimension(20) :: quark_qw
       real(8), parameter :: Pi=3.14159265D0, Nc=3.D0,nf=5.D0

       

       

      real(8) :: p0,p1,qhat,L,alpha_F,thetac,tmin,tform,kt,tauform
      real(8) :: fact,omega,omegas,abar,R,kts
      integer :: grid,q,n0

     

      COMMON /input/ p0,p1,qhat,L,alpha_F,thetac,abar

      COMMON /quenching/ quark_qw

      COMMON /radius/ R

     
      

      call gridvar(p,jac,jprim,split_F,split_A,nstep,p0,
     &     p1,2,nstep0,y0,num)

      fun_qw=quark_qw

   
      
      call  qw(fun_qw,y0,nstep0,4.5D0,abar,qhat,L)


c      fun_qw=0.8D0
      
      do k=1,nstep0
      
       do i=1,nstep

         do j=1,nstep

            q=jprim(i,j)
            
            tform=(2.*qhat/(p(j)*p(q)/p(i))**3/3.)**(0.25D0)

            tauform=2.D0/(p(j)*p(q)/p(i)*time0**2)
c            tauform=2.D0/(p(j)*time0**2)
c            tauform=0.D0
c            taud=1/(qhat*time0**2.D0/12.D0)**0.3333D0
c            tmin=max(thetac,tform)   ! coherence

            tmin=(1.D0/(L*p(j)*p(q)/p(i)/2.D0))**0.5D0 ! incoherence
c            tmin=(1.D0/(L*p(j)))**0.5D0 
            kt=time0*p(j)*p(q)/p(i)

            

              
            

            
            
              if (i<j) then
                 collision_F(i,j)=0.D0
                 collision_A(i,j)=0.D0
               
              elseif (time0 < tmin) then

                 collision_F(i,j)=0.D0
                 collision_A(i,j)=0.D0
                 
              elseif (kt< p0*1.01) then
                 
                 collision_F(i,j)=0.D0
                 collision_A(i,j)=0.D0
             
              else
                 omega=min(p(j),p(q))
                 omegas=qhat*L**2*(abar)**2.D0
c     &         *(1.D0-R**2/2.D0)
                 kts=sqrt(qhat*L)/R
c                 omegas=10.D0


                 if (omega < max(omegas,kts)) then
c                     fact=1.D0-exp(-R**2.D0)
                     fact=0.D0
                  else
                     fact=1.D0
                 endif

               collision_F(i,j)=
     &                (fun_A_in(j,k) * fun_F_in(q,k)* fun_qw(k)**2
     &          - fun_F_in(i,k))/time0
     &                /log(time0*p(j)*p(q)/p0/p(i))

c                  collision_F(i,j)= - fun_F_in(i,k)/time0
c     &                /log(time0*p(j)*p(q)/p0/p(i))

               
                  collision_A(i,j)=
     &              (fun_A_in(j,k) * fun_A_in(q,k)* fun_qw(k)**2
     &             *fact
     &          - fun_A_in(i,k)*fact)/time0
     &                /log(time0*p(j)*p(q)/p0/p(i))
               


               

            endif 
   
      
          enddo
      enddo      

      
      do i=1,nstep
         
         integral_F(i)=0.D0
         integral_A(i)=0.D0

*   gain term 

       if (i /= nstep) then
           do j=1,i

              integral_F(i)= integral_F(i)+ 0.5D0 *
     &         (jac(j)*split_F(i,j)*collision_F(i,j)
     &             +jac(j+1)*split_F(i,j+1)*collision_F(i,j+1))

                integral_A(i)= integral_A(i)+ 0.5D0 *
     &         (jac(j)*split_A(i,j)*collision_A(i,j)
     &         +jac(j+1)*split_A(i,j+1)*collision_A(i,j+1))
         
          enddo
  
          else
           do j=1,nstep-1 
              integral_F(nstep)= integral_F(nstep)+ 0.5D0 *
     &         (jac(j)*split_F(nstep,j)*collision_F(nstep,j)
     &             +jac(j+1)*split_F(nstep,j+1)*collision_F(nstep,j+1))
              
                integral_A(nstep)= integral_A(nstep)+ 0.5D0 *
     &         (jac(j)*split_A(nstep,j)*collision_A(nstep,j)
     &         +jac(j+1)*split_A(nstep,j+1)*collision_A(nstep,j+1))    
           enddo 
         endif


         fun_F_out(i,k)=integral_F(i)
         fun_A_out(i,k)=integral_A(i)
         
         
      enddo

      enddo


       

   
      
      end subroutine 


***************************************

*     grid and splitting function in pT definiton

***************************************
      
      
      subroutine gridvar(p_var,jaco,delt_ij,func_F,func_A,nstep,pmin,
     &  pmax,grid0,nstep0,p_var0,num)

   

      integer, intent(in) :: nstep,nstep0,grid0
      real(8), intent(in) :: pmin,pmax
      real(8), dimension(nstep), intent(out) :: p_var,jaco
      real(8), dimension(nstep0), intent(out) :: p_var0
      real(8), dimension(nstep,nstep), intent(out) :: func_F,func_A
      integer, dimension(nstep,nstep), intent(out) :: delt_ij
      integer, dimension(nstep0), intent(out) :: num
      
      integer  i,j
      real(8)  a,b,b0,delt_p,z
      real(8), dimension(nstep) :: y_var,dp_var
      real(8), dimension(nstep0) :: y_var0


      if(grid0 == 1) then        ! linear grid

         a =  pmax-pmin 
         b =  nstep * pmin - pmax 

         do i=1,nstep
            p_var(i)= ( a * i + b ) / (nstep-1)
            dp_var(i) = a / (nstep - 1)
            jaco(i)= abs(dp_var(i))
         enddo

         do i=1,nstep
            do j=1,nstep
               if(i > j) then
                  delt_ij(i,j)=i-j
               else
                  delt_ij(i,j)=0
                  
               endif
            enddo
         enddo   
         
      elseif(grid0 == 2) then    ! log grid :  y = log (1/x)
         
         a = - log (pmax/pmin) 
         b = - a * nstep
         b0= - a * nstep0
         
           
         do i=1,nstep
           y_var(i) = (a * i + b) / ( nstep - 1)
           p_var(i) = exp( - y_var(i) ) * pmax
           dp_var(i) = - (a / (nstep - 1) ) * p_var(i)

           jaco(i)= abs(dp_var(i))
      
        enddo
        
        do j=1,nstep0
           
           y_var0(j) = (a * j + b0) / ( nstep0 - 1)
           p_var0(j) = exp( - y_var0(j) ) * pmax
      
        enddo

        do j=1,nstep0
         num(j)= NINT( (nstep-1)*(a*j+b0)/(nstep0-1)/a-b/a)
        enddo
        
        
         do i=1,nstep
            do j=1,nstep
               delt_p = p_var(i)-p_var(j)
               if(delt_p >= pmin) then
                  delt_ij(i,j)=
     &          NINT(  (log(pmax/delt_p)*(nstep - 1)-b) / a )
               else
                  delt_ij(i,j)=0
                  
               endif
            enddo
         enddo   

      else
         write(*,*) 'invalid grid input, grid = 1,2,3'
      endif



** gluon-quark Splitting function 

       do i=1,nstep            
            do j=1,nstep
               z=p_var(j)/p_var(i)
               if(i > j) then
                
                 func_F(i,j) = (1.D0+(1.D0-z)**2)/p_var(j)
c                 func_F(i,j) = 2.D0/p_var(j)
                 func_A(i,j) = (1.D0-z*(1.D0-z))**2/p_var(j)/(1.D0-z)


               else
                  func_F(i,j) = 0.D0
                  func_A(i,j) = 0.D0

               endif                  
               
         enddo 
      enddo
      
      end subroutine


********************************************      
