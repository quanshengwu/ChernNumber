
  subroutine ham_qlayer2qlayer3(Hij)

     use para

     implicit none

! loop index
     integer :: iR

! index used to sign irvec     
     integer :: ia, ib

! new index used to sign irvec     
     integer :: new_ia, new_ib

! H00 Hamiltonian between nearest neighbour-quintuple-layers
! the factor 2 is induced by spin
     complex(Dp), intent(out) :: Hij(-ijmax:ijmax,-ijmax:ijmax,nband*2,nband*2)

     Hij= 0.0d0
     do iR=1, Nrpts
        ia= irvec(1, iR)
        ib= irvec(2, iR)

        ! 001
        new_ia=ia
        new_ib=ib
        
        if (abs(new_ib).le.ijmax.and.abs(new_ia).le.ijmax)then
           Hij(new_ia, new_ib , 1:Nband*2, 1:Nband*2 )&
           = HmnR(:,:,iR)
        endif

     enddo

     return
  end

