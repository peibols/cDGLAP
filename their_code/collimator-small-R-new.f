      program    collimator

      implicit none 

      integer, parameter :: n=1000, m=120,n0=20
      real(8), parameter :: Pi=3.14159265D0, Nc=3.D0,nf=5.D0
      real(8) :: CR,CF,CA,beta0,tscale,dtheta,tstep_F,tstep_A,time,taud
      real(8) :: qhat,L,R,omegac,thetac,p0,p1,alpha_F,alpha_A,abar,pmin
      real(8) :: omegas,jet_func,tag

      real(8), dimension(n)  :: p,jac
      real(8), dimension(n0)  :: y,f_F,f_A
      real(8), dimension(n,n0)  :: int_F_out,int_A_out,quark,gluon,
     & mint_F_out,mint_A_out,m_quark,m_gluon,
     &test_F_out,test_A_out,tmo_F_out,tmo_A_out,
     &test_in

      real(8), dimension(n0)  :: quark_qw,gluon_qw,pwr,pwr_F,pwr_A
      real(8), dimension(n,n) :: split_qq,split_gq,split_gg,split_qg

      
      integer, dimension(n,n) :: iprim
 
      integer, dimension(n0) :: num0

      integer :: grid,i,j,k,diag,initial
      
      character*(3)   arg
      character*(3)   arg2
      character*(3)   arg3


      COMMON /input/ p0,p1,qhat,L,alpha_F,thetac,abar

      COMMON /quenching/ quark_qw
      COMMON /radius/ R
      
      
      CF=(Nc**2-1.D0)/(2*Nc)
      CA=Nc

c  alpha_s(Q)=Pi/beta0/Log(Q/p0)

      beta0=(11.D0*Nc-2.D0*nf)/(6.D0)

c  IR cutoff chosen to match alpha_s(Mz)

      p0=0.09D0

      p1=5000.D0

      grid=2
     
      
c     prefactor of the collision integral
      

      tscale=1.D0/beta0


c     step in angular integral (time variable)
      
      dtheta=1.D-2
      

      tstep_F=dtheta*tscale

      tstep_A=dtheta*tscale

      
c     paramters

      qhat=1.5D0/5.D0
      L=17.5D0
      R=1.0D0


      alpha_F=CF/Pi*(Pi/beta0)

      alpha_A=CA/Pi*(Pi/beta0)

      abar=0.3
      
c     relevant scales

c     omega_c

      omegac= qhat * L**2 /12.D0
      omegas=qhat*L**2*(abar)**2.D0

c     theta_c

      thetac= sqrt(1.D0/omegac/L )
      pmin=(2.*qhat/3.)**0.333D0/R**(4.D0/3.D0)


       write(*,*) 'thetac, pmin, omegas, omegaR'
       write(*,*) thetac,pmin,omegas,sqrt(qhat*L)/R
***************************************

c      define grid for the integration over the splitting function
    
       call gridvar(p,jac,iprim,split_qq,split_gq,
     &split_gg,split_qg,n,p0,p1,2,n0,y,num0)


***************************************

*    Initial condition for dglap evolution at R=1
*   (from fits to Pythia)
* dn/dpt = A ( pt0 / pt ) ^ ( n + beta Log[pt/pt0]
* + gammaLog[pt/pt0]^2 )
* The values inside Params_... files correspond to the values of the 5 parameters, and are ordered as:
*  A
*  pt0
*  n
*  beta
*  gamma
*quark jet R=1
*0.22548980322510295
*48.55193282409811
*6.993359499440793
*-1.1854965558477715
*0.37816941197405374
*gluon jet R=1
*1.0580449529372458
*45.6161616210578
*7.475355563599871
*-1.146225037972374
*0.37486936745035404
**************
        initial=1

         if (initial==1) then

* quark jets
* power
      pwr_F= 6.993359499440793
     & - 1.1854965558477715 * Log(y/48.55193282409811)
     & + 0.37816941197405374 * Log(y/48.55193282409811)**2.
* factor
      f_F=0.22548980322510295
     &*(48.55193282409811/y)**(pwr_F)
* gluon jets
* power
      pwr_A=7.475355563599871
     & - 1.146225037972374 * Log(y/45.6161616210578)
     & + 0.37486936745035404 * Log(y/45.6161616210578)**2.
* factor
      f_A=1.0580449529372458
     &*(45.6161616210578/y)**(pwr_A)

       else if (initial==2) then
       
* quark jets
* power
              pwr_F =7.1199
     & - 1.51078 * Log(y/56.6076)
     & + 0.502422 * Log(y/56.6076)**2.
     
* factor
              f_F=0.000214679
     &*(56.6076/y)**(pwr_F)
* gluon jets
* power
              pwr_A=7.15593
     & - 1.25842 * Log(y/62.8043)
     & + 0.528596 * Log(y/62.8043)**2.


* factor
              f_A=0.000233974
     &*(62.8043/y)**(pwr_A)
          else
             write(*,*) 'incorrect initial condition'
          end if
************************


      write(*,*)  'testing the grid'
      write(*,*)  num0(1), split_qq(10,5),split_gg(10,5),num0(10)
      write(*,*)   y(5), p(num0(5))
     
c the power law and jet flavor
      
      pwr=pwr_A
      tag=0.D0   ! 1 for quark jet and 0for gluon jets
     
c calculate the quark and gluon quenching weights

      

      call qw(quark_qw,y,n0,pwr,alpha_F,qhat,L)
      call qw(gluon_qw,y,n0,pwr,alpha_A,qhat,L)


      write(*,*) 'test2:',y(15),quark_qw(15),pwr_F(15)


********************************************
*   initial condition
*******************************************
      
        do j=1,n0
           do i=1,n

            quark(i,j)=quark_qw(j)
            gluon(i,j)=gluon_qw(j)
            m_quark(i,j)=quark_qw(j)
            m_gluon(i,j)=gluon_qw(j)

           enddo 
        enddo
********************************************
*   angular evolution  (time)
*******************************************


      time=1.D-10



          do k=1,m

         

         if (mod(k,20) == 1) then
            write(arg,'(F3.1)')  time
            write(arg2,'(F3.1)')  R
         open(10,file='collimator/collimator_gluon_'
     &//arg//'_'//arg2//'.dat'
     &     ,status='unknown')
            write(arg,'(F3.1)')  time
         open(20,file='jet-spectrum/med_gluon_spect_'
     &//arg//'_'//arg2//'.dat'
     &     ,status='unknown')

              


         
         
         do j=1,n0

            diag=num0(j)

         jet_func=tag* f_F(j) * m_quark(diag,j)
     &   + (1.-tag)*f_A(j)* m_gluon(diag,j)


             write(10,*)    y(j), tag*quark(diag,j)
     &                    +(1.-tag)*gluon(diag,j)
             write(20,*)    y(j), jet_func
            


          enddo

            
        endif


          
          call coll(int_F_out,int_A_out,mint_F_out,
     &mint_A_out,quark,gluon,m_quark,m_gluon,n,n0,
     &time,pwr)

 

         do j=1,n0
           do i=1,n


              quark(i,j) = quark(i,j)+int_F_out(i,j)*tstep_F
              gluon(i,j) = gluon(i,j)+int_A_out(i,j)*tstep_A

             if (time<R) then
                  m_quark(i,j)=quark(i,j)
                  m_gluon(i,j)=gluon(i,j)
             else
                       
                
            m_quark(i,j) = m_quark(i,j)+mint_F_out(i,j)*tstep_F
            m_gluon(i,j) = m_gluon(i,j)+mint_A_out(i,j)*tstep_A


              endif
            
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
      real(8), dimension(nstep), intent(in) :: x,power
      real(8), dimension(nstep) :: integral
      real(8), dimension(nstep,mstep) :: spect,u
      real(8), dimension(mstep) :: w

      real(8), intent(in) :: alpha,length,qh
      integer, intent(in) :: nstep
     
      
      real(8) :: Pi,z,y,omegas,omegac,omegaBH,nu,R,temp,a,omega,ang

      COMMON /radius/ R

      Pi=3.14159265D0

      

      omegas=max(0.25**2*qh*length**2/1.5
     &,sqrt(qh*length/2.)/R)
      omegac=qh*length**2/2
      omegaBH=0.7D0
      z=omegaBH/omegas
      z=0.001D0
      y=omegac/omegas


      temp=0.7D0

      ang=1.-R**2*4./Pi**2



      
      do j=1,mstep
         
         w(j)=j*7.D0/mstep

         
         
      enddo

       do i=1,nstep
          do j=1,mstep
             omega= w(j)*x(i)/power(i)
             u(i,j)=sqrt(omegac*power(i)/(2*(x(i))*w(j)))
             
             if (u(i,j) > sqrt(omegac/2.D0/omegaBH)) then
                spect(i,j)=0.D0
             else if (omega > omegas) then
                spect(i,j)=0.D0
             else
c            spect(i,j) =log ( (cosh(u(i,j)))**2 - (sin(u(i,j)))**2 )
           spect(i,j) =log( (cosh(u(i,j))*cos(u(i,j)))**2.
     &     + (sinh(u(i,j))*sin(u(i,j)))**2.)
     &       /max(log((2.*qh*omega)**0.25D0/0.09D0),1.22)
c     &    *exp(-2.*(R*x(i)*w(j)/power(i))**2/qh/length)
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
     &     ( spect(i,j)* (1.D0 - exp(- w(j)*ang))/j
     &    + spect(i,j+1)* (1.D0 - exp(-w(j+1)*ang)) / (j+1) )
            

             

            
          enddo
        enddo

c        func=integral

       func =
     & exp(-2.*alpha*integral
     &      - 0.*qh/temp*power/x*length)
c       func =exp(- 0.3 * integral)
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
      subroutine coll(fun_F_out,fun_A_out,mom_F_out,mom_A_out,
     &fun_F_in,fun_A_in,mom_F_in,mom_A_in,nstep,nstep0,time0,
     &power)

      real(8), intent(in) :: time0
      real(8), dimension(nstep,nstep0), intent(in) ::
     &  fun_F_in,fun_A_in,mom_F_in,mom_A_in
      real(8), dimension(nstep,nstep0), intent(out) ::
     &  fun_F_out,fun_A_out,mom_F_out,mom_A_out
      real(8), dimension(nstep0), intent(in) ::
     &  power

      real(8), dimension(nstep) :: p,jac,gain
      real(8), dimension(nstep0) :: fun_qw,y0
      real(8), dimension(nstep,nstep) :: split_qq,split_gq,
     &split_gg,split_qg
      real(8), dimension(nstep,nstep) :: collision_A,collision_F,
     &mcoll_A,mcoll_F

       real(8), dimension(nstep) :: integral_A,integral_F,
     &mint_A,mint_F
       integer, dimension(nstep,nstep) :: jprim
       integer, dimension(nstep0) :: num

       
       real(8), parameter :: Pi=3.14159265D0, Nc=3.D0,nf=5.D0

       

       

      real(8) :: p0,p1,qhat,L,alpha_F,thetac,tmin,tform,kt,tauform
      real(8) :: fact,mfact,omega,omegas,abar,R,kts,zz,n_q,n_g
      integer :: grid,q,n0

     

      COMMON /input/ p0,p1,qhat,L,alpha_F,thetac,abar



      COMMON /radius/ R

     
       

      call gridvar(p,jac,jprim,split_qq,split_gq,split_gg,
     &split_qg,nstep,p0,p1,2,nstep0,y0,num)


      
      do k=1,nstep0
      
       do i=1,nstep

         do j=1,nstep
            

            q=jprim(i,j)

            zz=p(j)/p(i)
            
            tform=(2.*qhat/(p(j)*p(q)/p(i))**3/3.)**(0.25D0)

            tauform=2.D0/(p(j)*p(q)/p(i)*time0**2)

c            tmin=min(R,max(thetac,tform))  ! coherence
            tmin=max(thetac,tform)    ! coherence
c            tmin=(1.D0/(L*p(j)*p(q)/p(i)/2.D0))**0.5D0 ! incoherence

            kt=time0*p(j)*p(q)/p(i)

      
            
              if (i<j) then
                 collision_F(i,j)=0.D0
                 collision_A(i,j)=0.D0
                 mcoll_F(i,j)=0.D0
                 mcoll_A(i,j)=0.D0
               
              elseif (time0 < tmin) then

                 collision_F(i,j)=0.D0
                 collision_A(i,j)=0.D0
                 mcoll_F(i,j)=0.D0
                 mcoll_A(i,j)=0.D0
                 
              elseif (kt< p0*1.01) then
                 
                collision_F(i,j)=0.D0
                 collision_A(i,j)=0.D0
                 mcoll_F(i,j)=0.D0
                 mcoll_A(i,j)=0.D0
             
              else
                 omega=min(p(j),p(q))
                 omegas=qhat*L**2*(abar)**2.D0
c                 omegas=0.1
c     &         *(1.D0-R**2/2.D0)
                 kts=sqrt(qhat*L)/R




                 if (omega < max(omegas,kts)) then
c
                     fact=0.D0
                  else
                     fact=1.D0
                 endif

** Collision integral for the quenching weights

                collision_F(i,j)=split_qq(i,j)*fact
     &          *(fun_F_in(j,k) * fun_A_in(q,k)
     &          - fun_F_in(i,k))
     &    /time0/max(log(time0*p(j)*p(q)/p0/p(i)),1.22)
          

           
               
                  collision_A(i,j)=fact*(split_gg(i,j)
     &              *(fun_A_in(j,k) * fun_A_in(q,k)
     &             - fun_A_in(i,k))
     &           +split_qg(i,j)
     &              *(fun_F_in(j,k) * fun_F_in(q,k)
     &             - fun_A_in(i,k)))
     &    /time0/max(log(time0*p(j)*p(q)/p0/p(i)),1.22)
             
** Collision integral for the moments

     
      n_g=power(k)-1.


             mcoll_F(i,j)=fact*(split_qq(i,j)*
     &  (mom_F_in(j,k) *(fun_A_in(q,k)-1.)*zz**(n_g)
     &   +(mom_F_in(j,k)-mom_F_in(i,k))*zz**(n_g)
     &     - mom_F_in(i,k)*(1.-zz**(n_g)))+split_qq(i,j)*
     & ((fun_F_in(j,k)-1.) * mom_A_in(q,k)*(1.-zz)**(n_g)
     &    +mom_A_in(q,k)*(1.-zz)**(n_g))
     &    )/time0/max(log(time0*p(j)*p(q)/p0/p(i)),1.22)

               
            mcoll_A(i,j)=fact*(split_gg(i,j)*
     &  ( mom_A_in(j,k) * (fun_A_in(q,k)-1.)*zz**(n_g)
     & +(fun_A_in(j,k)-1.) *mom_A_in(q,k)*(1.-zz)**(n_g)
     &  +(mom_A_in(j,k)-mom_A_in(i,k))*zz**(n_g)
     &  + (mom_A_in(q,k)-mom_A_in(i,k))*(1.-zz)**(n_g)
     &  - mom_A_in(i,k)*(1.-zz**(n_g)-(1.-zz)**(n_g)))
     &  + split_qg(i,j)
     &   *(mom_F_in(j,k) * fun_F_in(q,k)*zz**(n_g)
     &  + fun_F_in(j,k) * mom_F_in(q,k)*(1.-zz)**(n_g)
     &         - mom_A_in(i,k))
     &   )/time0/max(log(time0*p(j)*p(q)/p0/p(i)),1.22)

            endif



          
          enddo
      enddo      


      do i=1,nstep
         
         integral_F(i)=0.D0
         integral_A(i)=0.D0
         mint_F(i)=0.D0
         mint_A(i)=0.D0
         
*   Integrating over z for the quenching weights and moments

       if (i /= nstep) then
           do j=1,i
               
              integral_F(i)= integral_F(i)+ 0.5D0 *
     &         (jac(j)*collision_F(i,j)
     &             +jac(j+1)*collision_F(i,j+1))
    
                integral_A(i)= integral_A(i)+ 0.5D0 *
     &         (jac(j)*collision_A(i,j)
     &         +jac(j+1)*collision_A(i,j+1))

         mint_F(i)= mint_F(i)+ 0.5D0 *
     &         (jac(j)*mcoll_F(i,j)
     &             +jac(j+1)*mcoll_F(i,j+1))

          mint_A(i)= mint_A(i)+ 0.5D0 *
     &         (jac(j)*mcoll_A(i,j)
     &         +jac(j+1)*mcoll_A(i,j+1))
         
          enddo
  
          else
           do j=1,nstep-1
               
              integral_F(nstep)= integral_F(nstep)+ 0.5D0 *
     &         (jac(j)*collision_F(nstep,j)
     &             +jac(j+1)*collision_F(nstep,j+1))
              
                integral_A(nstep)= integral_A(nstep)+ 0.5D0 *
     &         (jac(j)*collision_A(nstep,j)
     &         +jac(j+1)*collision_A(nstep,j+1))


          mint_F(nstep)= mint_F(nstep)+ 0.5D0 *
     &         (jac(j)*mcoll_F(nstep,j)
     &             +jac(j+1)*mcoll_F(nstep,j+1))
         
           mint_A(nstep)= mint_A(nstep)+ 0.5D0 *
     &         (jac(j)*mcoll_A(nstep,j)
     &         +jac(j+1)*mcoll_A(nstep,j+1))
           enddo 
         endif


         fun_F_out(i,k)=integral_F(i)
         fun_A_out(i,k)=integral_A(i)


         mom_F_out(i,k)=mint_F(i)
         mom_A_out(i,k)=mint_A(i)
         
         
      enddo

*   Integrating over z for the moments


*   enddo of the k = pT loop

      enddo


   
      
      end subroutine 


****************************************************

*     grid and splitting function definiton

*****************************************************
      
      
      subroutine gridvar(p_var,jaco,delt_ij,func_qq,func_gq,
     &func_gg,func_qg,nstep,pmin,pmax,
     &grid0,nstep0,p_var0,num)


* p_var : momentum
* jaco: jacobian for the p integrals
* delt_ij: finds the integer q corresponding to p(i)-p(j)
* nstep: size of the z (momentum) grid
* pmin: minimum momentum
* pmax: maximum momentum
* grid0: grid type, 0 for linear, 1 for logarithmic
* nstep0: size of the pT grid
* p_var0: pT
* num: returns p=pT
* power: power of the spectrum as function of pT

   
      real(8), parameter :: Pi=3.14159265D0, Nc=3.D0,nf=5.D0
      real(8)  CF,CA,TR
      integer, intent(in) :: nstep,nstep0,grid0
      real(8), intent(in) :: pmin,pmax
      real(8), dimension(nstep), intent(out) :: p_var,jaco
      real(8), dimension(nstep0), intent(out) :: p_var0
      real(8), dimension(nstep,nstep), intent(out) ::
     &func_qq,func_qg,func_gq,func_gg
      integer, dimension(nstep,nstep), intent(out) :: delt_ij
      integer, dimension(nstep0), intent(out) :: num
      
      integer  i,j
      real(8)  a,b,b0,delt_p,z
      real(8), dimension(nstep) :: y_var,dp_var
      real(8), dimension(nstep0) :: y_var0


      CA=Nc
      CF=(Nc**2.D0-1.D0)/Nc/2.D0
      
      TR=0.5D0

* Grid definition

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



** Splitting functions

       do i=1,nstep            
            do j=1,nstep
               z=p_var(j)/p_var(i)
               if(i > j) then
                
            func_gq(i,j) = CF*(1.D0+(1.D0-z)**2)/z/p_var(i)

            func_gg(i,j) =
     &           CA*(1.D0-z*(1.D0 - z))**2/(z*(1.D0-z))/p_var(i)
            func_qq(i,j)=CF*(1.D0 + z**2)/(1.D0-z)/p_var(i)
            func_qg(i,j)=TR*nf*((1.D0-z)**2+z**2)/p_var(i)


               else
            func_gq(i,j) = 0.D0

            func_gg(i,j) = 0.D0
                     
            func_qq(i,j) = 0.D0
            func_qg(i,j) = 0.D0

               endif                  
               
         enddo 
      enddo
      
      end subroutine


********************************************      
