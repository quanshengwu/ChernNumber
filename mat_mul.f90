! performs matrix-matrix multiply
! C=A*B
  subroutine mat_mul(m, A, B, C)
     
     use para, only : Dp, zzero, zone
     implicit none


     integer,intent(in) :: m

     complex(Dp) :: ALPHA
     complex(Dp) :: BETA 
 

     complex(Dp), intent(in)  :: A(m, m)
     complex(Dp), intent(in)  :: B(m, m)
     complex(Dp), intent(out) :: C(m, m)

     ALPHA= zone
     BETA = zzero

    !C(:,:)= zzero

     call ZGEMM('N','N',m,m,m,ALPHA, &
               &  A,m,B,m,BETA,C,m)

     return
  end subroutine mat_mul

! performs matrix-matrix multiply
! C=A*B'
  subroutine mat_mul_c(m, A, B, C)
     
     use para, only : Dp, zzero, zone
     implicit none


     integer,intent(in) :: m

     complex(Dp) :: ALPHA
     complex(Dp) :: BETA 
 

     complex(Dp), intent(in)  :: A(m, m)
     complex(Dp), intent(in)  :: B(m, m)
     complex(Dp), intent(out) :: C(m, m)

     ALPHA= zone
     BETA = zzero

    !C(:,:)= zzero

     call ZGEMM('N','C',m,m,m,ALPHA, &
               &  A,m,B,m,BETA,C,m)

     return
  end subroutine mat_mul_c

! performs matrix-matrix multiply
! C=A*B
! B is a diagnoal matrix
  subroutine mat_mul_diagB(m, A, B, C)
     
     use para, only : Dp
     implicit none


     integer,intent(in) :: m

     complex(Dp), intent(in)  :: A(m, m)
     complex(Dp), intent(in)  :: B(m)
     complex(Dp), intent(out) :: C(m, m)

     integer :: i
     integer :: j
 
     do i=1, m
        do j=1, m
           C(i,j)= A(i,j)*B(j)
        enddo
     enddo

     return
  end subroutine mat_mul_diagB

! performs matrix-matrix multiply
! C=A*B
! A is a diagnoal matrix
  subroutine mat_mul_diagA(m, A, B, C)
     
     use para, only : Dp
     implicit none


     integer,intent(in) :: m

     complex(Dp), intent(in)  :: A(m)
     complex(Dp), intent(in)  :: B(m, m)
     complex(Dp), intent(out) :: C(m, m)

     integer :: i
     integer :: j

     do i=1, m
        do j=1, m
           C(i,j)= A(i)*B(i,j)
        enddo
     enddo

     return
  end subroutine mat_mul_diagA


