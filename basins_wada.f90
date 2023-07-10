program basins_wada
  implicit none
  integer, parameter :: dp=kind(0.d0)
  integer :: i,j,k,condi1,condi2,Npts,saida,flag
  real(dp) :: pi,xk,yk,B_field, R_major,a_minor, M, L, v_parallel, alpha, gama, beta, w0,phi,E_r
  real(dp) :: x,y,x0,y0,aux1,aux2,xmin,xmax,ymin,ymax,q_factor,R1,R2,b, b2, g,qa,xe
  open(1,file="exit-1.dat")
  open(2,file="exit-2.dat")
  open(3,file="exit-3.dat")
  open(4,file="exits-grid.dat")
  !Declaração das variáveis
  pi=4.0d0*atan(1.0d0)
  Npts=int(1d5)
  condi1 = 1000 
  R_major = 10.d0/3.0d0; a_minor = 1.0d0; B_field = 1.0d0
  w0 = 16.364d0; M=15.d0; L=6.d0
  alpha = 3.570d0; beta = -7.921d0; gama = 4.132d0  
  phi= 9.84d-3 

  xmin = -0.5d0; xmax = 0.5d0
  ymin = 0.3d0; ymax = 1.0d0
  do j=1,condi1  
    x0 = xmin + dfloat(j-1)*(xmax-xmin)/dfloat(condi1-1)
    aux1 = x0     
    do i=1,condi1   
      y0 = ymin + dfloat(i-1)*(ymax-ymin)/dfloat(condi1-1)
      aux2 = y0  
      xk = x0
      yk = y0
      xe = 1.0d0/6.0d0 
      saida = 0
      k = 0
      !Iterate the map
      do while (k < Npts .and. saida == 0)
        y = yk + 4.0d0*pi*M*phi*dsin(xk*2.0d0*pi)/(w0)
        if(y<=1.0d0) then
	        q_factor = 5.0d0 - 6.3d0*y**2 + 6.3d0*y**3
	      else 
	        q_factor = 5.0*y
        end if
        v_parallel = -9.867d0 + 17.478d0*dtanh(10.1d0*y - 9.0d0)
        E_r = 3.0d0*alpha*y + 2.0d0*beta*dsqrt(dabs(y)) + gama
        R1 = (M-L*q_factor)*v_parallel/(w0*R_major*q_factor)
        R2 = -E_r*M/(w0*dsqrt(dabs(y)))
        x = xk + R1 + R2     
        x = dabs(x)
        x = x - int(x/(1.0d0)) 
        if (x .lt. -0.5d0) then
          x=x+1.d0
        endif
        if (x .gt. 0.5d0) then
          x=x-1.d0
        endif
        if (y>=1.0d0) then
          if (x>=-0.5d0 .and.x<-xe) then
            saida = 1
          endif
          if (x>=-xe.and.x<xe) then
            saida = 2
          endif
          if (x>=xe .and.x<=0.5d0) then
            saida = 3
          endif
        endif
        xk = x
        yk = y              
        k = k + 1
      end do  
    if (saida == 1) then
      write(1,"(f10.6,f10.6,f10.5,i10)") aux1, aux2,0.0d0, k
    endif  
    if (saida == 2) then
      write(2,"(f10.6,f10.6,f10.5,i10)") aux1, aux2, pi, k
    endif
    if (saida == 3) then
      write(3,"(f10.6,f10.6,f10.5,i10)") aux1, aux2, 2.0d0*pi, k
    endif
    write(4,"(f10.6,f10.6,i10)") aux1, aux2, saida
  end do    
  end do 
  close(1) 
  close(2)
  close(3)
  close(4)
  stop
end program basins_wada
