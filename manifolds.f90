real*8 f, g,xn, yn, xn1, yn1, pi, xi, yi,x0,y0,e1,e2,e3,v1,v2,v3,v4,alpha,beta,gamma,M,L
real*8 c1, x, y, xsela, ysela,xmin,xmax,ymin,ymax, vpara,w0,q,E,zz,col
integer*4 itera,condi1,condi2,i,j

  !Define the map 
  f(x,y)= (y + beta*dsin(2.0*pi*x))
  q(x,y) = 5.0d0 - 6.3d0*y**2 + 6.3d0*y**3
  E(x,y) = e1*y + e2*dsqrt(dabs(y)) + e3
  vpara(x,y) = v1 + v2*dtanh(v3*y + v4)
  g(x,y)= x + alpha*vpara(x,y)*(M/q(x,y) - L) + gamma*E(x,y)/dsqrt(dabs(y))
  !Define the area of interest to plot the manifolds, as well the grid.
  pi=4.0d0*atan(1.0d0)
  L=1000
  L=1000
  xmin = -0.5
  xmax = 0.5
  ymin = 0.8d0
  ymax = 1.0d0
  !Define the parameters of the map
  e1 = 10.7d0; e2 = -15.8d0; e3 = 4.13;
  M = 15.0d0; L = 6.0d0;
  v1 = -9.867d0; v2 = 17.47d0; v3 = 10.1d0; v4 = - 9.0d0
  alpha = 1.83d-2; gamma = - 9.16d-1
  beta = 0.10d0
  ! Define the number of iterations and if there is colission or not.
  col = 0.0d0
  itera = 10  

  open(1,file='stable.dat')
  open(2,file='unstable.dat')
  open(3,file='saddle.dat')

  do j=1,L     
    xi = xmin + dfloat(j-1)*(xmax-xmin)/dfloat(L-1) 
    do i=1,L
      yi = ymin + dfloat(i-1)*(ymax-ymin)/dfloat(condi2-1)  
      xn = xi
      yn = yi
      n = 1       
      do while ((n.lt.itera))
      CALL RANDOM_NUMBER(zz)    
      zz = zz - 0.5 ! Shift zz to [-0.5, 0.5]
      zz = zz*2.0d0*pi ! Scale phi to [-pi, pi]
      yn1=f(xn,yn) + col*dsin(zz)
      xn1=g(xn,yn1) + col*dcos(zz) 
      xn1 = dabs(xn1)
      xn1 = xn1 - int(xn1/(1.0d0))  
      if (xn1 .lt. -0.5d0) then
        xn1=xn1+1.d0
      endif
      if (xn1 .gt. 0.5d0) then
        xn1=xn1-1.d0
      endif

      xn=xn1
      yn=yn1  
      !We save the values for have iterations for plotting the saddle
      if(n.eq.itera/2)then
        xsela=xn
        ysela=yn
      endif
      
      n=n+1        
    enddo !Here is define the region where we calculate the manifolds.    
    if (yn .lt. 1.0d0 .and. yn .gt. 0.80d0)then
      if (yi .lt. 1.0d0 .and. yi .gt. 0.80d0)then
        if (ysela .lt. 1.0d0 .and. ysela .gt. 0.8d0)then
          if(n.eq.itera)then
            write(1,*) xi, yi, 2.0
          if (yn.lt.1.0d0 .and. yn.gt.0.80d0) then
            write(2,*) xn, yn, 2.0
          endif
        endif
          write(3,*) xsela, ysela
        endif
      endif
    endif
    enddo      
  enddo     
close(1)
close(2)
close(3)
end
