! calculate bulk's energy band using wannier TB method
  subroutine ek_bulk

     use para

     implicit none

     integer :: ik, i, j

     integer :: NN, knv3, nlines
     ! wave vector in 3d
     real(Dp) :: kstart(2), kend(2)

     real(dp) :: thetax
     real(dp) :: thetay
     real(Dp) :: W(ndim)

     real(dp) :: tt1, temp, kp(4,2)

     ! Hamiltonian of bulk system
     complex(Dp) :: Hamk_bulk(ndim,ndim) 

     real(dp), allocatable :: kpoint(:,:)
     
     ! eigen value of H
     real(dp), allocatable :: eigv(:,:)

     real(dp), allocatable :: k_len(:)

     complex(dp), allocatable :: dis (:,:)

     kp(1,:)=(/0.5d0, 0.5d0/)  ! R
     kp(2,:)=(/0.0d0, 0.0d0/)  ! Gamma
     kp(3,:)=(/0.5d0, 0.0d0/)  ! X
     kp(4,:)=(/0.5d0, 0.5d0/)  ! M
    
     nlines=3
     NN= 50
     knv3=NN*nlines
     allocate( kpoint(knv3, 2))
     allocate( eigv  (knv3, ndim))
     allocate( dis   (ndim, ndim))
     allocate( k_len (knv3))
     kpoint= 0d0
     eigv  = 0d0

     tt1=0d0
     k_len=0d0
     do j=1, nlines 
        do i=1, NN
           kstart= kp(j,:)
           kend  = kp(j+1,:)
           kpoint(i+(j-1)*NN,:)= kstart+ (kend-kstart)*dble(i-1)/dble(NN-1)
           
           temp= dsqrt((kp(j+1,1)- kp(j,1))**2 &
                 +(kp(j+1,2)- kp(j,2))**2)/dble(NN-1) 

           if (i.gt.1) then
              tt1=tt1+temp
           endif
           k_len(i+(j-1)*NN)= tt1
        enddo
     enddo


     do ik= 1, knv3
        print * , ik

        thetax = kpoint(ik, 1)*2d0*pi
        thetay = kpoint(ik, 2)*2d0*pi

        ! calculation bulk hamiltonian
        Hamk_bulk=0d0
        call ham_twist(thetax, thetay, Hamk_bulk)
        call disorder(dis)
        Hamk_bulk= Hamk_bulk+ dis

        W=0d0
        call eigensystem_c( 'N', 'U', ndim ,Hamk_bulk, W)

        eigv(ik, :)=W

     enddo

     open(unit=14, file='ektwist.dat')

    !do ik=1, knv3
    !   write(14, '(100000f19.9)')k_len(ik),eigv(ik, :)
    !enddo

     do i=1, Ndim
        do ik=1, knv3
           write(14, '(10f20.10)')k_len(ik), eigv(ik,i)
        enddo
        write(14, *)' '
     enddo

     do ik=1, Nlines
        write(14, *)k_len((ik-1)*NN+1), -5d0
        write(14, *)k_len((ik-1)*NN+1),  5d0
        write(14, *)' '
     enddo
     write(14, *)k_len(knv3), -5d0
     write(14, *)k_len(knv3),  5d0
     write(14, *)' '

     close(14)

   return
   end subroutine ek_bulk
