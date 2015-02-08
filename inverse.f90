  subroutine inv(ndim,Amat)

     implicit none

     integer,parameter :: dp=8

     integer           :: i
     integer           :: info

     integer,intent(in):: ndim

!    IPIV   : INTEGER. Array, DIMENSION at least max(1, n). The pivot indices that define
!    the permutation matrix P; row i of the matrix was interchanged with
!    row ipiv(i). Corresponds to the single precision factorization (if info=
!    0 and iter ≥ 0) or the double precision factorization (if info= 0 and
!    iter < 0).
     integer,allocatable      :: ipiv(:)

     complex(dp),parameter        :: zzero=(0.0d0,0.0d0)
     complex(dp),parameter        :: zone=(1.0d0,0.0d0)
     complex(dp),parameter        :: zi=(0.0d0,1.0d0)

!    Amat  :
!    Overwritten by the factors L and U from the factorization of A = P*L*U;
!    the unit diagonal elements of L are not stored.
!    If iterative refinement has been successfully used (info= 0 and iter≥
!    0), then A is unchanged.
!    If double precision factorization has been used (info= 0 and iter <
!    0), then the array A contains the factors L and U from the factorization
!    A = P*L*U; the unit diagonal elements of L are not stored.
     complex(dp),intent(inout):: Amat(ndim,ndim)

!    Bmat  :
!    Identity matrix
!    Overwritten by the solution matrix X for dgesv, sgesv,zgesv,zgesv.
     complex(dp),allocatable :: Bmat(:, :)

     allocate(ipiv(ndim))
     allocate(Bmat(ndim, ndim))

     ipiv= 0
     Bmat= zzero
     do i=1, ndim
        Bmat(i, i)= zone
     enddo

     call zgesv(ndim,ndim,Amat,ndim,ipiv,Bmat,ndim,info)

     if (info.ne.0) then
        print *,'something wrong with cpu zgesv','info=',info
        stop
     endif

     Amat=Bmat
     
     deallocate(ipiv)
     
     return
  end subroutine inv 

