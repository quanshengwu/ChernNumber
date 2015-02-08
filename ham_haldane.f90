! This subroutine is used to caculate Hamiltonian for 
! slab system with twist boundary
! with zigzag edge haldane model
! 

! History  
!        Dec/11th/2012 by Quansheng Wu
!        Dec/22th/2012 by Quansheng Wu

  subroutine ham_haldane(thetax, thetay, hamk_twist)
  
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
     integer :: i

     complex(Dp) :: expiphi

     complex(Dp) :: ratiox
     complex(Dp) :: ratioy

     complex(Dp), allocatable :: H00x(:, :)
     complex(Dp), allocatable :: H01x(:, :)

     allocate(H00x(Ny*2, Ny*2))
     allocate(H01x(Ny*2, Ny*2))
     H00x= zzero
     H01x= zzero

     !> construct H00x
     do i=1, Ny
        H00x((i-1)*2+1, (i-1)*2+1)= M
        H00x((i-1)*2+2, (i-1)*2+2)=-M
     enddo
     do i=1, Ny*2-1
        H00x(i, i+1)= t0
        H00x(i+1, i)= t0
     enddo
     expiphi= exp(zi*phi)


     do i=1, Ny/2-1
        H00x((i-1)*4+1, (i-1)*4+3)= t1*expiphi
        H00x((i-1)*4+2, (i-1)*4+4)= t1*expiphi
        H00x((i-1)*4+5, (i-1)*4+3)= t1*expiphi
        H00x((i-1)*4+6, (i-1)*4+4)= t1*expiphi
        H00x((i-1)*4+3, (i-1)*4+5)= t1*conjg(expiphi)
        H00x((i-1)*4+4, (i-1)*4+6)= t1*conjg(expiphi)
        H00x((i-1)*4+3, (i-1)*4+1)= t1*conjg(expiphi)
        H00x((i-1)*4+4, (i-1)*4+2)= t1*conjg(expiphi)
     enddo
     i= Ny/2
     H00x((i-1)*4+1, (i-1)*4+3)= t1*expiphi
     H00x((i-1)*4+2, (i-1)*4+4)= t1*expiphi
     H00x((i-1)*4+3, (i-1)*4+1)= t1*conjg(expiphi)
     H00x((i-1)*4+4, (i-1)*4+2)= t1*conjg(expiphi)

     ratioy= cos(thetay)+ zi* sin(thetay)
     !> twist boundary along y direction
     H00x((i-1)*4+4, 1)= H00x((i-1)*4+4, 1)+ t0* ratioy
     H00x((i-1)*4+4, 2)= H00x((i-1)*4+4, 2)+ t1*conjg(expiphi)* ratioy
     H00x((i-1)*4+3, 1)= H00x((i-1)*4+3, 1)+ t1*conjg(expiphi)* ratioy
     H00x(1, (i-1)*4+4)= H00x(1, (i-1)*4+4)+ t0* conjg(ratioy)
     H00x(2, (i-1)*4+4)= H00x(2, (i-1)*4+4)+ t1*expiphi* conjg(ratioy)
     H00x(1, (i-1)*4+3)= H00x(1, (i-1)*4+3)+ t1*expiphi* conjg(ratioy)

     !> construct H01x
     do i=1, Ny*2-1, 4
        H01x(i, i)= t1*conjg(expiphi)
        H01x(i+1, i)= t0
        H01x(i+1, i+1)= t1*expiphi
        H01x(i+1, i+3)= t1*conjg(expiphi)
        H01x(i+2, i+3)= t0
        H01x(i+2, i  )= t1*expiphi
        H01x(i+2, i+2)= t1*conjg(expiphi)
        if(i+4.le.2*Ny) H01x(i+2, i+4)= t1*expiphi
     enddo

     do i=1, Ny/2
        H01x(i*4, i*4)= t1*expiphi
     enddo
     i= Ny/2

     !> twist boundary along y direction
     H01x((i-1)*4+3, 1)= H01x((i-1)*4+3, 1)+ t1*expiphi*ratioy
     H01x(2, (i-1)*4+4)= H01x(2, (i-1)*4+4)+ t1*conjg(expiphi)*conjg(ratioy)

     hamk_twist=0.0d0 
     do i= 1, Nx
        hamk_twist((i-1)*2*Ny+1:i*2*Ny, (i-1)*2*Ny+1:i*2*Ny)=  H00x
     enddo

     do i= 1, Nx-1
        hamk_twist((i-1)*2*Ny+1:i*2*Ny, i*2*Ny+1:(i+1)*2*Ny)= H01x
        hamk_twist(i*2*Ny+1:(i+1)*2*Ny, (i-1)*2*Ny+1:i*2*Ny)= conjg(transpose(H01x))
     enddo

     !>>  twist boundary along x direction
     ratiox= cos(thetax)+ zi* sin(thetax)

     i1= Nx
     i2= 1 
     hamk_twist((i1-1)*2*Ny+1:i1*2*Ny, (i2-1)*2*Ny+1:i2*2*Ny)= &
     hamk_twist((i1-1)*2*Ny+1:i1*2*Ny, (i2-1)*2*Ny+1:i2*2*Ny)+  H01x*ratiox
     hamk_twist((i2-1)*2*Ny+1:i2*2*Ny, (i1-1)*2*Ny+1:i1*2*Ny)= &
     hamk_twist((i2-1)*2*Ny+1:i2*2*Ny, (i1-1)*2*Ny+1:i1*2*Ny)+ conjg(transpose(H01x)*ratiox)

! check hermitcity

     do i1=1, Ndim
     do i2=i1+1, Ndim
        if(abs(hamk_twist(i1,i2)-conjg(hamk_twist(i2,i1))).ge.1e-6)then
          write(*,*)'there are something wrong with hamk_twist'
          stop
        endif 
     enddo
     enddo

     return
  end


  subroutine ham_haldane2(thetax, thetay, hamk_twist)
  
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
     integer :: i

     complex(Dp) :: expiphi

     complex(Dp), allocatable :: H00x(:, :)
     complex(Dp), allocatable :: H01x(:, :)

     allocate(H00x(Ny*2, Ny*2))
     allocate(H01x(Ny*2, Ny*2))
     H00x= zzero
     H01x= zzero

     !> construct H00x
     do i=1, Ny
        H00x((i-1)*2+1, (i-1)*2+1)= M
        H00x((i-1)*2+2, (i-1)*2+2)=-M
     enddo
     do i=1, Ny*2-1
        H00x(i, i+1)= t0
        H00x(i+1, i)= t0
     enddo
     expiphi= exp(zi*phi)


     ! next nearest neighbour
     do i=1, Ny/2-1
        H00x((i-1)*4+1, (i-1)*4+3)= t1*expiphi
        H00x((i-1)*4+2, (i-1)*4+4)= t1*expiphi
        H00x((i-1)*4+5, (i-1)*4+3)= t1*expiphi
        H00x((i-1)*4+6, (i-1)*4+4)= t1*expiphi
        H00x((i-1)*4+3, (i-1)*4+5)= t1*conjg(expiphi)
        H00x((i-1)*4+4, (i-1)*4+6)= t1*conjg(expiphi)
        H00x((i-1)*4+3, (i-1)*4+1)= t1*conjg(expiphi)
        H00x((i-1)*4+4, (i-1)*4+2)= t1*conjg(expiphi)
     enddo
     i= Ny/2
     H00x((i-1)*4+1, (i-1)*4+3)= t1*expiphi
     H00x((i-1)*4+2, (i-1)*4+4)= t1*expiphi
     H00x((i-1)*4+3, (i-1)*4+1)= t1*conjg(expiphi)
     H00x((i-1)*4+4, (i-1)*4+2)= t1*conjg(expiphi)

     !> construct H01x
     do i=1, Ny*2-1, 4
        H01x(i, i)= t1*conjg(expiphi)
        H01x(i+1, i)= t0
        H01x(i+1, i+1)= t1*expiphi
        H01x(i+1, i+3)= t1*conjg(expiphi)
        H01x(i+2, i  )= t1*expiphi
        H01x(i+2, i+2)= t1*conjg(expiphi)
        H01x(i+2, i+3)= t0
        if(i+4.le.2*Ny) H01x(i+2, i+4)= t1*expiphi
     enddo

     do i=1, Ny/2
        H01x(i*4, i*4)= t1*expiphi
     enddo

     hamk_twist=0.0d0 
     do i= 1, Nx
        hamk_twist((i-1)*2*Ny+1:i*2*Ny, (i-1)*2*Ny+1:i*2*Ny)= H00x
     enddo

     do i= 1, Nx-1
        hamk_twist((i-1)*2*Ny+1:i*2*Ny, i*2*Ny+1:(i+1)*2*Ny)= H01x
        hamk_twist(i*2*Ny+1:(i+1)*2*Ny, (i-1)*2*Ny+1:i*2*Ny)= conjg(transpose(H01x))
     enddo

! check hermitcity

     do i1=1, Ndim
     do i2=i1+1, Ndim
        if(abs(hamk_twist(i1,i2)-conjg(hamk_twist(i2,i1))).ge.1e-6)then
          write(*,*)'there are something wrong with hamk_twist'
          stop
        endif 
     enddo
     enddo

     return
  end

  subroutine ham_xy( X, Y)
  
     use para
     implicit none

     ! Hamiltonian of slab system
     complex(Dp),intent(out) :: X(ndim, ndim) 
     complex(Dp),intent(out) :: Y(ndim, ndim) 

     ! loop index  
     integer :: i

     complex(Dp), allocatable :: H00x(:, :)
     complex(Dp), allocatable :: H00y(:, :)

     allocate(H00x(Ny*2, Ny*2))
     allocate(H00y(Ny*2, Ny*2))
     H00x= zzero
     H00y= zzero

     !> construct H00x
     do i=1, Ny
        H00x((i-1)*2+1, (i-1)*2+1)= 1d0
        H00x((i-1)*2+2, (i-1)*2+2)= 1d0
        H00y((i-1)*2+1, (i-1)*2+1)= dble(i)
        H00y((i-1)*2+2, (i-1)*2+2)= dble(i)
     enddo
    

     X=0.0d0 
     Y=0.0d0 
     do i= 1, Nx
        X((i-1)*2*Ny+1:i*2*Ny, (i-1)*2*Ny+1:i*2*Ny)= dble(i)*H00x
        Y((i-1)*2*Ny+1:i*2*Ny, (i-1)*2*Ny+1:i*2*Ny)= H00y
     enddo

     return
  end
