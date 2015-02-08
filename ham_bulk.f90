! This subroutine is used to caculate Hamiltonian for 
! slab system with open  boundary

! History  
!        Dec/11th/2012 by Quansheng Wu

  subroutine ham_bulk(thetax, thetay, hamk_bulk)
  
     use para
     implicit none

     ! wave vector in 2d
     real(Dp), intent(in) :: thetax
     real(Dp), intent(in) :: thetay

     ! Hamiltonian of slab system
     complex(Dp),intent(out) ::hamk_bulk(ndim, ndim) 

     ! loop index  
     integer :: i1
     integer :: i2
     integer :: j1
     integer :: j2

     complex(Dp) :: ratiox
     complex(Dp) :: ratioy

     complex(Dp), allocatable :: Hij(:, :, :, :)

     allocate(Hij(-ijmax:ijmax, -ijmax:ijmax, nband*2, nband*2))

     call ham_qlayer2qlayer3(Hij)


     hamk_bulk=0.0d0 
     do i1= 1, Nx
     do j1= 1, Ny
     do i2= 1, Nx   
     do j2= 1, Ny   
       if (abs(i2-i1).le.ijmax .and. abs(j2-j1).le.ijmax ) then
           hamk_bulk( (i1-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                      (i1-1)*Nband*2*Ny+j1*Nband*2,       &
                      (i2-1)*Nband*2*Ny+(j2-1)*Nband*2+1: &
                      (i2-1)*Nband*2*Ny+j2*Nband*2 )      &
           = Hij(i2-i1, j2-j1, 1:nband*2, 1:nband*2)

       endif 
     enddo
     enddo
     enddo
     enddo

     return

     !>>  bulk boundary
     ratiox= cos(thetax)+ zi* sin(thetax)
     ratioy= cos(thetay)+ zi* sin(thetay)

     j1= Ny
     j2= 1 
     do i1= 1, Nx
     do i2= 1, Nx   

        if (abs(i2-i1).le.ijmax)then
           hamk_bulk( (i1-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                      (i1-1)*Nband*2*Ny+j1*Nband*2,       &
                      (i2-1)*Nband*2*Ny+(j2-1)*Nband*2+1: &
                      (i2-1)*Nband*2*Ny+j2*Nband*2 )      &
           = Hij(i2-i1, 1, 1:nband*2, 1:nband*2)*ratioy
        endif

     enddo
     enddo

     j1= 1 
     j2= Ny
     do i1= 1, Nx
     do i2= 1, Nx   

        if (abs(i2-i1).le.ijmax)then
           hamk_bulk( (i1-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                      (i1-1)*Nband*2*Ny+j1*Nband*2,       &
                      (i2-1)*Nband*2*Ny+(j2-1)*Nband*2+1: &
                      (i2-1)*Nband*2*Ny+j2*Nband*2 )      &
           = Hij(i2-i1, -1, 1:nband*2, 1:nband*2)*conjg(ratioy)
        endif

     enddo
     enddo

 

     i1= Nx
     i2= 1 
     do j1= 1, Ny
     do j2= 1, Ny   

        if (abs(j2-j1).le.ijmax)then
           hamk_bulk( (i1-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                      (i1-1)*Nband*2*Ny+j1*Nband*2,       &
                      (i2-1)*Nband*2*Ny+(j2-1)*Nband*2+1: &
                      (i2-1)*Nband*2*Ny+j2*Nband*2 )      &
           = Hij(1, j2-j1, 1:nband*2, 1:nband*2)*ratiox
        endif

     enddo
     enddo

     i1= 1 
     i2= Nx
     do j1= 1, Ny
     do j2= 1, Ny   

        if (abs(j2-j1).le.ijmax)then
           hamk_bulk( (i1-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                      (i1-1)*Nband*2*Ny+j1*Nband*2,       &
                      (i2-1)*Nband*2*Ny+(j2-1)*Nband*2+1: &
                      (i2-1)*Nband*2*Ny+j2*Nband*2 )      &
           = Hij(-1, j2-j1, 1:nband*2, 1:nband*2)*conjg(ratiox)
        endif

     enddo
     enddo
 
! check hermitcity

     do i1=1, Ndim
     do i2=1, Ndim
        if(abs(hamk_bulk(i1,i2)-conjg(hamk_bulk(i2,i1))).ge.1e-6)then
          write(*,*)'there are something wrong with hamk_bulk'
          stop
        endif 
     enddo
     enddo

  return
  end
