! This subroutine is used to caculate Hamiltonian for 
! slab system with twist boundary

! History  
!        Dec/11th/2012 by Quansheng Wu

  subroutine ham_twist(thetax, thetay, hamk_twist)
  
     use para
     implicit none

     ! wave vector in 2d
     real(Dp), intent(in) :: thetax
     real(Dp), intent(in) :: thetay

     ! Hamiltonian of slab system
     complex(Dp),intent(out) ::hamk_twist(ndim, ndim) 

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


     hamk_twist=0.0d0 
     do i1= 1, Nx
     do j1= 1, Ny
     do i2= 1, Nx   
     do j2= 1, Ny   
       if (abs(i2-i1).le.ijmax .and. abs(j2-j1).le.ijmax ) then
           hamk_twist((i1-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                      (i1-1)*Nband*2*Ny+j1*Nband*2,       &
                      (i2-1)*Nband*2*Ny+(j2-1)*Nband*2+1: &
                      (i2-1)*Nband*2*Ny+j2*Nband*2 )      &
           = Hij(i2-i1, j2-j1, 1:nband*2, 1:nband*2)

       endif 
     enddo
     enddo
     enddo
     enddo

!    return
     !>>  twist boundary
     ratiox= cos(thetax)+ zi* sin(thetax)
     ratioy= cos(thetay)+ zi* sin(thetay)

     i1= Nx
     i2= 1 
     do j1= 1, Ny
     do j2= 1, Ny   

        if (abs(j2-j1).le.ijmax)then
           hamk_twist((i1-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                      (i1-1)*Nband*2*Ny+j1*Nband*2,       &
                      (i2-1)*Nband*2*Ny+(j2-1)*Nband*2+1: &
                      (i2-1)*Nband*2*Ny+j2*Nband*2 )      &
           = Hij(1, j2-j1, 1:nband*2, 1:nband*2)*ratiox+  &
           hamk_twist((i1-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                      (i1-1)*Nband*2*Ny+j1*Nband*2,       &
                      (i2-1)*Nband*2*Ny+(j2-1)*Nband*2+1: &
                      (i2-1)*Nband*2*Ny+j2*Nband*2 )      
        endif

     enddo
     enddo

     i1= 1 
     i2= Nx
     do j1= 1, Ny
     do j2= 1, Ny   

        if (abs(j2-j1).le.ijmax)then
           hamk_twist((i1-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                      (i1-1)*Nband*2*Ny+j1*Nband*2,       &
                      (i2-1)*Nband*2*Ny+(j2-1)*Nband*2+1: &
                      (i2-1)*Nband*2*Ny+j2*Nband*2 )      &
           = Hij(-1, j2-j1, 1:nband*2, 1:nband*2)*conjg(ratiox) + &
           hamk_twist((i1-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                      (i1-1)*Nband*2*Ny+j1*Nband*2,       &
                      (i2-1)*Nband*2*Ny+(j2-1)*Nband*2+1: &
                      (i2-1)*Nband*2*Ny+j2*Nband*2 )      
        endif

     enddo
     enddo

     j1= Ny
     j2= 1 
     do i1= 1, Nx
     do i2= 1, Nx   

        if (abs(i2-i1).le.ijmax)then
           hamk_twist((i1-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                      (i1-1)*Nband*2*Ny+j1*Nband*2,       &
                      (i2-1)*Nband*2*Ny+(j2-1)*Nband*2+1: &
                      (i2-1)*Nband*2*Ny+j2*Nband*2 )      &
           = Hij(i2-i1, 1, 1:nband*2, 1:nband*2)*ratioy+  &
           hamk_twist((i1-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                      (i1-1)*Nband*2*Ny+j1*Nband*2,       &
                      (i2-1)*Nband*2*Ny+(j2-1)*Nband*2+1: &
                      (i2-1)*Nband*2*Ny+j2*Nband*2 )      
        endif

     enddo
     enddo

     j1= 1 
     j2= Ny
     do i1= 1, Nx
     do i2= 1, Nx   

        if (abs(i2-i1).le.ijmax)then
           hamk_twist((i1-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                      (i1-1)*Nband*2*Ny+j1*Nband*2,       &
                      (i2-1)*Nband*2*Ny+(j2-1)*Nband*2+1: &
                      (i2-1)*Nband*2*Ny+j2*Nband*2 )      &
           = Hij(i2-i1, -1, 1:nband*2, 1:nband*2)*conjg(ratioy)+ &
           hamk_twist((i1-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                      (i1-1)*Nband*2*Ny+j1*Nband*2,       &
                      (i2-1)*Nband*2*Ny+(j2-1)*Nband*2+1: &
                      (i2-1)*Nband*2*Ny+j2*Nband*2 )      
        endif

     enddo
     enddo

 
 
! check hermitcity

     do i1=1, Ndim
     do i2=1, Ndim
        if(abs(hamk_twist(i1,i2)-conjg(hamk_twist(i2,i1))).ge.1e-6)then
          write(*,*)'there are something wrong with hamk_twist'
          stop
        endif 
     enddo
     enddo

  return
  end
