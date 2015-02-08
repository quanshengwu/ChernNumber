! complex version on cpu
! a subroutine to calculate eigenvector and eigenvalue
  subroutine eigensystem_c(JOBZ,UPLO,N,A,W)

     use para, only : Dp
     implicit none

!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
     character*1, intent(in) :: JOBZ

!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
     character*1, intent(in) :: UPLO

!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
     integer,   intent(in) :: N

!  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
!          On entry, the Hermitian matrix A.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of A contains the
!          upper triangular part of the matrix A.  If UPLO = 'L',
!          the leading N-by-N lower triangular part of A contains
!          the lower triangular part of the matrix A.
!          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!          orthonormal eigenvectors of the matrix A.
!          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!          or the upper triangle (if UPLO='U') of A, including the
!          diagonal, is destroyed.

     complex(Dp),intent(inout) :: A(N,N)

!  W       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, the eigenvalues in ascending order.

     real(Dp), intent(inout) :: W(N)
    
     integer :: info

     integer :: lwork

     real(Dp),allocatable ::  rwork(:)

     complex(Dp),allocatable :: work(:)

     lwork= 16*N
     allocate(rwork(lwork))
     allocate( work(lwork))


     W= 0.0d0
     info= 0
     call zheev( JOBZ, UPLO, N, A, N,  &
              W, work, lwork, rwork, info )

     return
  end subroutine


! complex version on cpu
! a subroutine to calculate eigenvector and eigenvalue
  subroutine eigensystem_zge(JOBZ, N, A, W)

     use para, only : Dp
     implicit none

!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
     character*1, intent(in) :: JOBZ

!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
     integer,   intent(in) :: N

!  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
!          On entry, the Hermitian matrix A.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of A contains the
!          upper triangular part of the matrix A.  If UPLO = 'L',
!          the leading N-by-N lower triangular part of A contains
!          the lower triangular part of the matrix A.
!          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!          orthonormal eigenvectors of the matrix A.
!          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!          or the upper triangle (if UPLO='U') of A, including the
!          diagonal, is destroyed.

     complex(Dp),intent(inout) :: A(N,N)

!  W       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, the eigenvalues in ascending order.

     complex(Dp), intent(inout) :: W(N)
    
     integer :: info

     integer :: lwork

     real(Dp),allocatable ::  rwork(:)

     complex(Dp),allocatable :: work(:)
     complex(Dp),allocatable :: VL(:, :)
     complex(Dp),allocatable :: VR(:, :)

     lwork= 16*N
     allocate(rwork(lwork))
     allocate( work(lwork))
     allocate( VL(N, N))
     allocate( VR(N, N))


     info=0
     W=0.0d0
     call zgeev(JOBZ, 'N', N, A, N,  &
              W, VL, N, VR, N, work, lwork, rwork, info )
     A= VR

     return
  end subroutine


