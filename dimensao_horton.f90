program dimensao
  
  implicit none
  integer, parameter :: dp=kind(0.0d0)
  real(dp) :: s,p,pi,x0,x,xmin,xmax,y0,y,ymin,ymax,aux1,aux2
  real(dp) :: theta,f,epsi,raio,x1,y1,x2,y2,x3,y3,r,xe, q1, q2,q3
  real(dp) :: q_factor,E_r,v_parallel,e1,e2,e3,v1,v2,v3,v4,alpha,beta,gamma,M,U
  integer :: i,j,k,l,n,nci,nstep,condi,saida,saida1
  integer :: flag1, flag2
  logical :: incert
  open(1,file="50-5.dat")
  
  !Declaração das variáveis
  pi = 4.0d0*atan(1.0d0)
  nci = 0
  f = 0.0d0
  condi = 1d6
  nstep = 2d6
  epsi = 1.0d-5

  e1 = 10.7d0; e2 = -15.8d0; e3 = 4.13;
  M = 15.0d0; U = 6.0d0;
  v1 = -9.867d0; v2 = 17.48d0; v3 = 10.1d0; v4 = - 9.0d0
  alpha = 1.83d-2; gamma = - 9.16d-1
  beta = 4.0*M*pi*7.65d-3/16.364!0.056686928d0


  !Intervalos de x e y
  xmin = -0.5d0
  xmax = 0.5d0 
  ymin = 0.8d0
  ymax = 1.0d0 

  call random_seed()
  
  do n=1,10
     
     !Condições iniciais
     i = 0
     nci = 0 !para calcular o erro de f
     
     do while (i < condi)
        
        !cond. inicias aleatórias
        call random_number(r)
        x0 = xmin + (xmax - xmin)*r
        aux1 = x0
        y0 = ymin + (ymax - ymin)*r
        aux2 = y0
        
        flag1 = 0
        saida = 0
        j = 0
        
        do while (j < nstep .and. saida == 0)
           j = j + 1
            y1 = y0 + beta*dsin(2.0d0*pi*x0)
            E_r = e1*y1 + e2*dsqrt(dabs(y1)) + e3
            if(y1<=1.0d0) then
	            q_factor = 5.0d0 - 6.3d0*y1**2 + 6.3d0*y1**3
	            else 
	                q_factor = 5.0*y1
            end if
            v_parallel = v1 + v2*dtanh(v3*y1 + v4)
            x1 = dmod(x0 + alpha*v_parallel*(M/q_factor - U) + gamma*E_r/dsqrt(dabs(y1))+0.5d0, 1.0d0)-0.5d0
                
            if (y1>1.0d0) then
               if (x1>=-0.5d0 .and. x1<0.0d0) then
                  saida = 1
               endif
               if (x1>=0.0d0 .and. x1<=0.5d0) then
                  saida = 2
               endif
            endif

           y0 = y1
           x0 = x1
        enddo
        
        if (saida == 1) then
              flag1 = 1
          else if (saida == 2) then
              flag1 = 1
          else
              flag1 = 0
        endif
        
        incert = .false.
        
        !     y2 = aux2 - epsi
        if (saida == 1) then
           
           do k=1,2
              
              call random_number(theta)
              theta = theta*2*pi
              call random_number(raio)
              raio = epsi*raio
              
              y2 = aux2 + raio*cos(theta)
              x2 = aux1 + raio*sin(theta)
              
              flag2 = 0
              saida1 = 0
              l = 0
              
              do while (l < nstep .and. saida1 == 0)   !loop de iteracao do mapa
                 l = l + 1
                 
                  y3 = y2 + beta*dsin(2.0d0*pi*x2)
                  E_r = e1*y3 + e2*dsqrt(dabs(y3)) + e3
                  if(y3<=1.0d0) then
	                  q_factor = 5.0d0 - 6.3d0*y3**2 + 6.3d0*y3**3
	               else 
	                  q_factor = 5.0*y3
                  end if
                  v_parallel = v1 + v2*dtanh(v3*y3 + v4)
                  x3 = dmod(x2 + alpha*v_parallel*(M/q_factor - U) + gamma*E_r/dsqrt(dabs(y3))+0.5d0, 1.0d0)-0.5d0
                
                  if (y3>1.0d0) then
                     if (x3>=-0.5d0 .and. x3<0.0d0) then
                        saida1 = 1
                     endif
                  if (x3>=0.0d0 .and. x3<=0.5d0) then
                     saida1 = 2
                  endif
               endif
                 y2 = y3
                 x2 = x3
                 
              enddo
              

              if (saida1 == 1) then
                   flag2 = 1
                else if (saida1 == 2) then
                   flag2 = 1
                else 
                   flag2 = 0
              endif
              
              if((flag1 .ne. flag2) .and. (incert .eqv. .false.)) incert = .true.
              
           enddo
        endif
        
        if(incert .eqv. .true.) nci = nci + 1
        
        if (saida == 1) i = i + 1
        
     enddo
     
     f = (dfloat(nci)/dfloat(condi))
     
     write(1,"(2(f20.15))") epsi,f
     
  enddo
  
  close(1)
  stop
end program dimensao

                