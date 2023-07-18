! Calculate the escape basins for a grid of L x L uniformly
! distributed initial conditions for the horton map with parameters from TCABr
! with N iterations.
! Author: Leonardo Costa de Souza
! Last modified: 30/06/2023
! Compile: gfortran stdmap.f90 bacia-horton.f90 -o exe.x  or use ifort (faster in general)
program basins_horton
   implicit none
   integer, parameter :: dp=kind(0.0d0)
   integer :: i,j,k,condi1,condi2,Npts,saida,flag
   real(dp) :: pi,xk,yk,e1,e2,e3,v1,v2,v3,v4,alpha,beta,gamma,M,L
   real(dp) :: x,y,x0,y0,aux1,aux2,xmin,xmax,ymin,ymax,q_factor,E_r,v_parallel
   open(1,file="exit-1.dat")
   open(2,file="exit-2.dat")
   open(3,file="entropy-1.dat")

   !Parameters
   pi=4.0d0*atan(1.0d0)
   N=int(1d5) ! Max number of iterations
   L = 1000; ! Grid
   !Parameters of the machine normalized
   e1 = 10.7d0; e2 = -15.8d0; e3 = 4.13;
   M = 15.0d0; L = 6.0d0;
   v1 = -9.867d0; v2 = 17.48d0; v3 = 10.1d0; v4 = - 9.0d0
   alpha = 1.83d-2; gamma = - 9.16d-1
   beta = 0.10d0
   ! Region of interest in the phase space
   xmin = -0.5; xmax = 0.5
   ymin = 0.375d0; ymax = 1.0d0
 
   do j=1,L     
      x0 = xmin + dfloat(j-1)*(xmax-xmin)/dfloat(L - 1)
      aux1 = x0 !save the initial condtion of x
      do i=1,L 
         y0 = ymin + dfloat(i-1)*(ymax-ymin)/dfloat(condi2-1)
         aux2 = y0 !save the initial condition of  y
         xk = x0
         yk = y0       
         saida = 0
         k = 0
         !Iterate the map testing if the particle has escape
         do while (k < N .and. saida == 0)
            y = yk + beta*dsin(2.0d0*pi*xk)
            E_r = e1*y + e2*dsqrt(dabs(y)) + e3
            if(y<=1.0d0) then
	            q_factor = 5.0d0 - 6.3d0*y**2 + 6.3d0*y**3
	            else 
	                q_factor = 5.0*y
            end if
            v_parallel = v1 + v2*dtanh(v3*y + v4)
            x = dmod(xk + alpha*v_parallel*(M/q_factor - L) + gamma*E_r/dsqrt(dabs(y))+0.5d0, 1.0d0)-0.5d0
            ! Now we test if the particle has reach one of the exits 
            if (y>1.0d0) then 
               if (x>=-0.5d0 .and. x<0.0d0) then
                  saida = 1
               endif
               if (x>=0.0d0 .and. x<=0.5d0) then
                  saida = 2
               endif
            endif
            xk = x
            yk = y              
            k = k + 1
         end do  
         ! Here, we write the values of x and y (initial conditions) if they have escape in 1 if the exit is 1 
         ! or 2 if it is the second exit
         ! Combine the two files and with them we plot the basins and also the escape time.
         if (saida == 1) then
            write(1,"(f10.6,f10.6,f10.5,i10)") aux1, aux2, pi, k
         endif  
         if (saida == 2) then
            write(2,"(f10.6,f10.6,f10.5,i10)") aux1, aux2, 2.0d0*pi, k
         endif
         write(3,"(f10.6,f10.6,i10)") aux1, aux2, saida
      end do    
   end do
close(1) 
close(2)
close(3)
stop
end program basin_horton
