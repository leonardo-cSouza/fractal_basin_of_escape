! Calculate if a given system is Wada or not using the Grdi-approach
! In this exemple the system is the horton map with parameters from TCABr
! with N iterations.
! Author: Leonardo Costa de Souza
! Last modified: 30/06/2023
program grid_wada
  implicit none
  integer, parameter :: dp=kind(0.d0)
  real(dp) :: x(1000,1000), y(1000,1000), pi, G0, G1, G2, G3, xi
  real(dp) :: xf, yi, yf, xtest, ytest, HISTO, dist, s, p
  integer(dp) :: i, j ,k
  integer :: M(1000,1000), sinal0, sinal1, sinal2, sinal3, i2, j2
  integer :: q, n, nmax, verifica, saida, wada, W(1000,1000)
  integer :: B(1000,1000), F(1000,1000), qmax, qmin, pegou, npts
  common/parameters/ pi, s, p, Npts
  
  pi=dacos(-1.d0)
  s = 2.0d0*pi
  p = 2.5d0
  Npts  = 1000000
  
  qmin=1
  qmax=20
	
!       File where is read the 1000x1000 initial conditions of the basins
  open(1,file='exits.dat')
!       Output file of W2 and W3 as a function of the refinement q
  open(5,file='W2W3_mt.dat')
!	      Output file hit the histogram of the reclassification of points G2 and G3
  open(6,file='HISTOGRAM_mt.dat')
	
  do i=1, 1000
     do j=1, 1000
        W(i,j)=0
        B(i,j)=0
        F(i,j)=0
     enddo
  enddo
  
	
!******  Set the matrices of positions and values of exits from the basins file
  do k=1, 1000000
     i=1+int(dfloat(k-1)/1000.d0)
     j=k-1000*(i-1)
     read(1,*) x(i,j), y(i,j), M(i,j)
  enddo
  G0=0.d0
  !*************************************************
  !******  Variation of the step q **********************
  do q=qmin, qmax
     HISTO=0.D0
     G1=0.d0
     G2=0.d0
     G3=0.d0  
     !       Sweeping the elements of the Basin Matrix
     do i=2, 999
        do j=2, 999
           
           if(W(i,j).eq.1) G3=G3+1.d0
           if(B(i,j).eq.1) G1=G1+1.d0
           !print*,G3
           
           if((F(i,j).eq.1).OR.(q.eq.qmin))then  
              sinal0=0
              sinal1=0
              sinal2=0
              sinal3=0
              pegou=0
              !       Compare the point M(I,J) with his 8 close neighbours 
              do i2=(i-1), (i+1)   
                 do j2=(j-1), (j+1)
                    if(M(i2,j2).eq.0) sinal0=1
                    if(M(i2,j2).eq.1) sinal1=1
                    if(M(i2,j2).eq.2) sinal2=1
                    if(M(i2,j2).eq.3) sinal3=1
                    
                    !       If the central point belong to a different basin of one of his neighbours,
                    !       we save the position of this two points to do a line
                    
                    dist=dsqrt((dfloat(i)-dfloat(i2))**2.d0+(dfloat(j)-dfloat(j2))**2.d0) 
                    if((M(i2,j2).ne.M(i,j)).and.(pegou.eq.0))then
                       xi=x(i,j)
                       yi=y(i,j)
                       xf=x(i2,j2)
                       yf=y(i2,j2)
                       if(dist.lt.1.4d0) pegou=1
                    endif
                 enddo
              enddo
              !       End of the sweepin of the point M(I,J)
              
              !       If the point is internal of a basin
              if(((sinal1+sinal2+sinal3).eq.1).and.(sinal0.eq.0))then
                 G1=G1+1.d0
                 B(i,j)=1
              endif
              
              !       if it is in the bondary of three basins
              if(((sinal1+sinal2+sinal3).eq.3).and.(sinal0.eq.0))then
                 G3=G3+1.d0
                 W(i,j)=1
              endif
              
              !       If it is in the boundary of two basins
              if(((sinal1+sinal2+sinal3).eq.2).and.(sinal0.eq.0))then
                 
                 !       In the step Q=1, this point belong to the group G2
                 if(q.eq.1) then
                    G2=G2+1.d0
                    F(i,j)=1
                 endif
                 
                 !*****************In the steps Q>1 There are refinaments that can reclassified as G3
                 if((q.gt.1).and.(W(i,j).eq.0)) then
                    verifica=0
                    
                    !       A given "verifica", that correspont to the value of exits necessary to 
                    !       this point be classified as G3
                    
                    if(sinal1.eq.0) verifica=1
                    if(sinal2.eq.0) verifica=2
                    if(sinal3.eq.0) verifica=3
                    
                    nmax=int((2.d0**(dfloat(q)-1.d0))/2.d0)
                    n=1
                    wada=0
                    do while ((n.le.nmax).and.(wada.eq.0))
                       
                       !       subroutine that determine the point (X,Y) to be test in this				   
                       call stepq (xi,yi,xf,yf,q,n,xtest,ytest)
                       saida=0
                       !       
                       !        Subroutine the determine the exit for the point (Xtest,Ytest)
                       call Map(xtest,ytest, saida)
                       
                       !       If the exit of this point (Xtest,Ytest) is the same value as "VERIFICA"
                       !       the point G2 becomes G3
                       if(saida.eq.verifica) wada=1
                       n=n+1
                    enddo
                    
                    if(wada.eq.1)then
                       G3=G3+1.d0
                       W(i,j)=1
                       F(i,j)=0
                       HISTO=HISTO+1.D0
                    else
                       if(saida.ne.0)then
                          G2=G2+1.d0
                          F(i,j)=1
                       else
                          F(i,j)=0
                          G0=G0+1.d0
                       endif
                    endif
                 endif
                 
                 !***********************************************************************************************************
                 !*****************End of the if of the points in the boundary of two basins 
              endif
           endif
        enddo
        
     enddo
!*       End of the sweeping of the elements of the basin Matrix
     
     write(5,*)  Q, G2/(G2+G3), G3/(G2+G3)
     WRITE(6,*)  Q, HISTO
     
     !*       File where is save the points G1
     open(2,file='out1.dat',status='REPLACE')
     write(2,*) 'refinamento', q
     !*       file where is save the points G2
     open(3,file='out2.dat',status='REPLACE')
     write(3,*) 'refinamento', q
     !*       File where is save the points G3
     open(4,file='out3.dat',status='REPLACE')
     write(4,*) 'refinamento', q
     do i=1, 1000
        do j=1, 1000
           if(B(i,j).eq.1) write(2,*) i, j
           if(F(i,j).eq.1) write(3,*) i, j
           if(W(i,j).eq.1) write(4,*) i, j
        enddo
     enddo
     close(2)
     close(3)
     close(4)
     
  enddo
  !******  End of the variation of the step q
  
  stop 
end program grid_wada
	
subroutine stepq (xi, yi, xf, yf, q, n, x, y)
  integer, parameter :: dp=kind(0.d0)
  real(dp) :: xi, yi, xf, yf, x, y
  integer :: q, n
!  real*8 xi, yi, xf, yf, x, y
!  integer*4 q, n

  
  if(xi.ne.xf)then
     x=xi+(2.d0*dfloat(n)-1.d0)/(2.d0**(dfloat(q)-1.d0))*(xf-xi)
     y=yi+(yf-yi)*(x-xi)/(xf-xi)
  else
     x=xi
     y=yi+(2.d0*dfloat(n)-1.d0)/(2.d0**(dfloat(q)-1.d0))*(yf-yi)
  endif
end subroutine stepq

Subroutine Map (xk,yk,saida)
  integer, parameter :: dp=kind(0.d0)
  real(dp) :: pi, s, p, xk, yk, x, y, xe
  integer :: Npts, saida, k, sinal
 ! real*8 pi, s, p, xk, yk, x, y, xe
 ! integer*4 Npts, saida, k, sinal
  common/parameters/ pi, s, p, Npts
  k =0
  sinal=0
  xe=-1.d0
  do while ((k.lt.Npts).and.(sinal.eq.0))
     x = xk + s*yk - p*dexp(-yk)*dcos(xk) + s*(dlog(dabs(dcos(xk - p*dexp(-yk)*dcos(xk)))) - dlog(dabs(dcos(xk))))
     y = yk + dlog(dabs(dcos(xk - p*dexp(-yk)*dcos(xk)))) - dlog(dabs(dcos(xk))) 
     
     if (y.lt. 0.0) then
        xe = xk - (x - xk)*yk/(y - yk)
        sinal=1
     endif
     
     do while (x .lt. 0.0d0)
        x = x + 2.0d0*pi
     enddo
     
     do while (x .gt. 2.0d0*pi)
        x = x - 2.0d0*pi
     enddo
     
     xk = x
     yk = y
     k=k+1
  enddo
  
  do while (xe.lt.0.0)
     xe = xe + 2.0*pi
  enddo
  
  do while (xe.gt.2.0*pi)
     xe = xe - 2.0*pi
  enddo
   if (sinal == 1) then
      if (xe>=0.d0 .and. xe<2.d0*pi/3.d0) then
         saida = 1
      endif    
      if(xe>=2.d0*pi/3.d0 .and. xe<4.d0*pi/3.d0) then
         saida = 2
      endif
      if(xe>4.d0*pi/3.d0 .and. xe<=2.0d0*pi) then
         saida = 3
         endif
      else 
         saida = 0
      endif
  if(saida.gt.3) saida=3
end Subroutine Map
