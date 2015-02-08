!>>> this method comes from 2010 EPL 92 67004 
! Disordered topological insulators via C*-algebras
! http://iopscience.iop.org/0295-5075/92/6/67004)
! T. A. Loring and M. B. Hastings
! this code is written by QuanSheng Wu on 12th Dec 2012
   subroutine chern_realspace(Ef, res)

      use para
      implicit none

      !>> inout parameter
      real(dp), intent(in) :: Ef
      real(dp), intent(out) :: res


      !>> local variables
      integer :: i
      integer :: j
      integer :: k
      integer :: it
      integer :: nfill

      real(dp) :: thetax
      real(dp) :: thetay

      real(dp), allocatable :: Eig(:)
      real(dp), allocatable :: xi(:)
      real(dp), allocatable :: yi(:)


      complex(dp), allocatable :: logy(:)
      complex(dp), allocatable :: ham(:, :)
      complex(dp), allocatable :: dis(:, :)
      complex(dp), allocatable :: PX(:, :)
      complex(dp), allocatable :: PY(:, :)
      complex(dp), allocatable :: TT(:, :)
      complex(dp), allocatable :: P (:, :)
      complex(dp), allocatable :: zeig(:)

      allocate(ham(ndim, ndim))
      allocate(dis(ndim, ndim))
      allocate(PX(ndim, ndim))
      allocate(PY(ndim, ndim))
      allocate(P (ndim, ndim))
      allocate(TT(ndim, ndim))
      allocate(eig(ndim))
      allocate(zeig(ndim))
      allocate(logy(ndim))
      allocate(xi(Ndim))
      allocate(yi(Ndim))
      PX= zzero
      PY= zzero
      P = zzero

      it=0
      do i=1, Nx
         do j=1, Ny
            do k=1, Nband*2
               it= it+1
               Xi(it)= dble(i)-1d0
               Yi(it)= dble(j)-1d0
            enddo
         enddo
      enddo

      Xi= Xi/dble(Nx)*2d0*pi
      Yi= Yi/dble(Ny)*2d0*pi

      ! PX= exp(i 2 pi X)
      ! PY= exp(i 2 pi Y)
      do i=1, ndim
         PX(i, i)= cos(xi(i))+ zi* sin(xi(i))
         PY(i, i)= cos(yi(i))+ zi* sin(yi(i))
      enddo
      
      thetax = 0d0
      thetay = 0d0

      ! prepare hamiltonian
     !call ham_bulk(thetax, thetay, ham)
      call ham_haldane(thetax, thetay, ham)
      call disorder(dis)
      ham= ham+ dis


      call eigensystem_c('V', 'U', Ndim, ham, eig)

      nfill= 0
      do i=1, ndim
         if (eig(i).le.Ef) nfill= nfill+1
      enddo

      call zgemm('N', 'C', Ndim, Ndim, nfill, zone, ham, &
                  Ndim, ham, Ndim, zzero, P, Ndim)


      ! P exp(i X 2 pi) P
      call mat_mul(Ndim, P, PX, TT)
      call mat_mul(Ndim, TT, P, PX)

      ! P exp(i Y 2 pi) P
      call mat_mul(Ndim, P, PY, TT)
      call mat_mul(Ndim, TT, P, PY)

      ! PX*PY*PX'*PY'
      call mat_mul  (Ndim, PX,  PY, TT)
      call mat_mul_c(Ndim, TT,  PX, ham)
      call mat_mul_c(Ndim, ham, PY, TT)


      call eigensystem_zge('N', Ndim, TT, zeig)


      res= 0d0
      do i=1, nfill
         res= res+ aimag(log(zeig(i)))/2d0/pi
      enddo

      return

   end subroutine chern_realspace

!>>> this method comes from  Ref: 2011 J. Phys. A: Math. Theor. 44 113001
! Disordered topological insulators: a non-commutative geometry perspective
! Emil Prodan
! this code is written by QuanSheng Wu on 17th Dec 2012
   subroutine chern_prodan(Ef, res)

      use para
      implicit none

      !>> inout parameter
      real(dp), intent(in) :: Ef
      real(dp), intent(out) :: res


      !>> local variables
      integer :: i
      integer :: j
      integer :: k
      integer :: it
      integer :: nfill

      ! this parameter comes from Eq (120, 121, 122) in Ref
      integer :: Q

      real(dp) :: thetax
      real(dp) :: thetay

      ! Delta= 2*pi/Nx
      real(dp) :: Delta

      real(dp) :: rt1

      real(dp), allocatable :: Eig(:)
      real(dp), allocatable :: xi(:)
      real(dp), allocatable :: yi(:)

      ! BB(1)=1/(2*Delta)
      complex(dp), allocatable :: BB(:)

      ! CC=Aij^-1*BB
      complex(dp), allocatable :: CC(:)

      ! Aij= j^{2i-1}
      complex(dp), allocatable :: Aij(:, :) 

      complex(dp), allocatable :: ham(:, :)
      complex(dp), allocatable :: dis(:, :)
      complex(dp), allocatable :: PX(:, :)
      complex(dp), allocatable :: PY(:, :)
      complex(dp), allocatable :: expx(:)
      complex(dp), allocatable :: TT1(:, :)
      complex(dp), allocatable :: TT2(:, :)
      complex(dp), allocatable :: xppx(:, :)
      complex(dp), allocatable :: yppy(:, :)
      complex(dp), allocatable :: P (:, :)

      ! the expansion order
      if (Nx.ge.16) then
         Q= 8 
      else
         Q= Nx/2
      endif

      allocate(ham(ndim, ndim))
      allocate(dis(ndim, ndim))
      allocate(PX(ndim, Q))
      allocate(PY(ndim, Q))
      allocate(P (ndim, ndim))
      allocate(BB (Q))
      allocate(CC (Q))
      allocate(Aij(Q, Q))
      allocate(TT1(ndim, ndim))
      allocate(TT2(ndim, ndim))
      allocate(xppx(ndim, ndim))
      allocate(yppy(ndim, ndim))
      allocate(eig(ndim))
      allocate(xi(Ndim))
      allocate(yi(Ndim))
      allocate(expx(ndim))
      PX= zzero
      PY= zzero
      P = zzero

      it=0
      do i=1, Nx
         do j=1, Ny
            do k=1, Nband*2
               it= it+1
               Xi(it)= dble(i)
               Yi(it)= dble(j)
            enddo
         enddo
      enddo


      ! PX= exp(I *i * X*Delta)
      ! PY= exp(I *i * Y*Delta)
      Delta= 2d0*pi/dble(Nx)
      do i=1, Q
         do j=1, Ndim
            rt1= dble(i)*Xi(j)*Delta
            PX(j, i)= cos(rt1)+ zi* sin(rt1)
            rt1= dble(i)*Yi(j)*Delta
            PY(j, i)= cos(rt1)+ zi* sin(rt1)
         enddo
      enddo

      ! > Aij
      do i=1, Q
         do j=1, Q
            Aij(i,j)= dble(j)**(2*i-1) 
         enddo
      enddo
      BB= zzero
      BB(1)= 1d0/2d0/Delta

      ! > calculate CC= Aij^-1*BB
      call inv(Q, Aij)
      i=1
      j=1
      call ZGEMV('N',Q,Q,zone,Aij,Q,BB,i,zzero,CC,j)
      
      thetax = 0d0
      thetay = 0d0

      ! prepare hamiltonian
     !call ham_bulk(thetax, thetay, ham)
     !call ham_twist(thetax, thetay, ham)
      call ham_haldane(thetax, thetay, ham)
      call disorder(dis)
      ham= ham+ dis

      call eigensystem_c('V', 'U', Ndim, ham, eig)

      nfill= 0
      do i=1, ndim
         if (eig(i).le.Ef) nfill= nfill+1
      enddo

      call zgemm('N', 'C', Ndim, Ndim, nfill, zone, ham, &
                  Ndim, ham, Ndim, zzero, P, Ndim)


      ! > Eq.(130)
      xppx= zzero
      do i=1, Q
         expx= conjg(PX(:,i))
         call mat_mul_diagA(Ndim, expx, P, TT1)
         expx= PX(:,i)
         call mat_mul_diagB(Ndim, TT1, expx, TT2)
         xppx= xppx+ zi* CC(i)* TT2

         expx= PX(:,i)
         call mat_mul_diagA(Ndim, expx, P, TT1)
         expx= conjg(PX(:,i))
         call mat_mul_diagB(Ndim, TT1, expx, TT2)
         xppx= xppx- zi* CC(i)* TT2
      enddo


      yppy= zzero
      do i=1, Q
         expx= conjg(PY(:,i))
         call mat_mul_diagA(Ndim, expx, P, TT1)
         expx= PY(:,i)
         call mat_mul_diagB(Ndim, TT1, expx, TT2)
         yppy= yppy+ zi* CC(i)* TT2

         expx= PY(:,i)
         call mat_mul_diagA(Ndim, expx, P, TT1)
         expx= conjg(PY(:,i))
         call mat_mul_diagB(Ndim, TT1, expx, TT2)
         yppy= yppy- zi* CC(i)* TT2
      enddo

      !> P[xppx, yppy], Eq.(129)
      call mat_mul(Ndim, xppx, yppy, TT1)
      call mat_mul(Ndim, yppy, xppx, TT2)
      TT1= TT1- TT2
      call mat_mul(Ndim, P, TT1, TT2)


      res= 0d0
      do i=1, Ndim
         res= res+ aimag(TT2(i, i))*2d0*pi/dble(Nx)/dble(Ny)
      enddo

      return

   end subroutine chern_prodan

!>>> this method comes from  Ref: 
! Journal of the Physical Society of Japan
! Vol. 74, No. 6, June, 2005, pp. 1674–1677
! Chern Numbers in Discretized Brillouin Zone:
! Efficient Method of Computing (Spin) Hall Conductances
! Takahiro FUKUI, Yasuhiro HATSUGAI and Hiroshi SUZUKI
! This code is written by QuanSheng Wu on 18th Dec 2012
! chern number in Discretized Brillouin Zone
!  subroutine chern_DBZ(Ef, res)

!     use para
!     implicit none


!     integer :: ik, i, j
!
!     integer :: NN, knv3

!     ! number of occupied states
!     integer :: nfill
!
!     real(dp) :: thetax
!     real(dp) :: thetay

!     integer, allocatable :: fill(:)
!     real(Dp), allocatable :: W(:)

!     ! Hamiltonian of bulk system
!     complex(Dp), allocatable :: Ham(:, :) 
!     complex(dp), allocatable :: dis (:,:)
!
!     integer , allocatable :: ikpoint(:,:)
!     real(dp), allocatable :: kpoint(:,:)
!     
!     complex(dp), allocatable :: zx(:)
!     complex(dp), allocatable :: zx(:)

!     complex(dp), allocatable :: U1_mat(:, :)
!     complex(dp), allocatable :: U2_mat(:, :)
!     complex(dp), allocatable :: detU1(:, :)
!     complex(dp), allocatable :: detU2(:, :)
!     complex(dp), allocatable :: eigvec(:, :, :)


!     NN= 10 
!     knv3= NN*NN
!     allocate( kpoint(knv3, 2))
!     allocate(ikpoint(knv3, 2))
!     allocate(W(ndim))
!     allocate(ham(ndim, ndim))
!     allocate(dis(ndim, ndim))
!     allocate(zx(ndim))
!     allocate(zy(ndim))
!     allocate(U1_mat(ndim, ndim))
!     allocate(U2_mat(ndim, ndim))
!     allocate(eigvec(ndim, ndim, knv3))
!     allocate(detU1(NN, NN))
!     allocate(detU2(NN, NN))
!     allocate(fill(knv3))
!     
!     ! > asign kpoint in [0, 1)*[0, 1)*4*pi*pi
!     ik= 0
!     do i=1, NN
!        do j=1, NN
!           ik= ik+ 1
!           kpoint(ik, 1)= dble(i-1)*/dble(NN)
!           kpoint(ik, 2)= dble(j-1)*/dble(NN)
!           ikpoint(ik, 1)= i
!           ikpoint(ik, 2)= j
!        enddo
!     enddo

!     kpoint= kpoint*2d0*pi
!
!     call disorder(dis)

!     ! >> prepare eigenvectors
!     do ik= 1, knv3
!
!        thetax = kpoint(ik, 1)
!        thetay = kpoint(ik, 2)
!
!        ! calculation bulk hamiltonian
!        Ham=0d0
!        call ham_twist(thetax, thetay, Ham)
!        Ham= Ham+ dis
!
!        W=0d0
!        call eigensystem_c( 'N', 'U', ndim ,Ham, W)

!        eigvec(:, :, ik)= ham

!        nfill= 0
!        do i=1, ndim
!           if (W(i).le.Ef)nfill= nfill+ 1
!        enddo
!        fill(ik)= nfill

!     enddo

!     ! >> calculate U1(k)
!     do ik= 1, knv3
!        nfill= fill(ik)
!        do i= 1, nfill
!           do j= 1, nfill
!              zx= ham(:, i)
!              zy= ham(:, j)
!              U_mat(i, j)= ZDOTC(ndim, ZX, INCX, ZY, INCY)                
!           enddo 
!        enddo 

!        ! det{U}
!        call det(ndim, nfill, U_mat, det1) 
!        i= ikpoint(ik, 1)
!        j= ikpoint(ik, 1)
!        detU(i, j)= det1/abs(det1)
!     enddo


!     return

!  end subroutine chern_DBZ

!>>> this method comes from  
! Ref: 
!   1. PRB 59 8144(1999)
!   2. JOURNAL OF APPLIED PHYSICS 112, 044306 (2012)
! construct on 18th Dec 2012
   subroutine chern_kubo(Ef, res)

      use para
      implicit none

      !>> inout variables
      real(dp), intent(in) :: Ef
      real(dp), intent(out) :: res

      integer :: ik, i, j
 
      integer :: NN, knv3

      ! number of occupied states
      integer :: nfill
 
      real(dp) :: thetax
      real(dp) :: thetay

      real(dp) :: En
      real(dp) :: Em

      complex(dp) :: z1
      complex(dp) :: z2
      complex(dp) :: z3

      real(Dp), allocatable :: W(:)
      real(dp), allocatable :: sigma(:)

      ! Hamiltonian of bulk system
      complex(Dp), allocatable :: Ham(:, :) 
      complex(dp), allocatable :: dis (:,:)
 
      real(dp), allocatable :: kpoint(:,:)
      
      complex(dp), allocatable :: zn(:)
      complex(dp), allocatable :: zm(:)
      complex(dp), allocatable :: vec1(:)

      complex(dp), allocatable :: vx(:, :)
      complex(dp), allocatable :: vy(:, :)
      complex(dp), allocatable :: vx0(:, :)
      complex(dp), allocatable :: vy0(:, :)

      complex(dp), allocatable :: H01x(:, :)
      complex(dp), allocatable :: H01y(:, :)

      ! vector multiply vector v1*v2
      complex(dp), external :: zdotc

      allocate(H01x(nband*2, nband*2))
      allocate(H01y(nband*2, nband*2))


      NN= 10
      knv3= NN*NN
      allocate( kpoint(knv3, 2))
      allocate(W(ndim))
      allocate(ham(ndim, ndim))
      allocate(dis(ndim, ndim))
      allocate(zn(ndim))
      allocate(zm(ndim))
      allocate(vec1(ndim))
      allocate(vx(ndim, ndim))
      allocate(vy(ndim, ndim))
      allocate(vx0(ndim, ndim))
      allocate(vy0(ndim, ndim))
      allocate(sigma(ndim))
      
      ! > asign kpoint in [0, 1)*[0, 1)*4*pi*pi
      ik= 0
      do i=1, NN
         do j=1, NN
            ik= ik+ 1
            kpoint(ik, 1)= dble(i-1)/dble(NN)
            kpoint(ik, 2)= dble(j-1)/dble(NN)
         enddo
      enddo

      kpoint= kpoint*2d0*pi
 
      call disorder(dis)

      ! > calculate H01x, H01y
      do ik=1, nrpts
         i = irvec(1, ik)
         j = irvec(2, ik)
         if (i .eq. 1 .and. j .eq. 0)H01x= Hmnr(:, :, ik)
         if (i .eq. 0 .and. j .eq. 1)H01y= Hmnr(:, :, ik)
      enddo
      call velocity2(ndim, vx0, vy0)

      ! >> prepare eigenvectors
      do ik= 1, knv3
 
         thetax = kpoint(ik, 1)
         thetay = kpoint(ik, 2)
 
         ! calculation bulk hamiltonian
         Ham=0d0
         call ham_twist(thetax, thetay, Ham)
         Ham= Ham+ dis
 
         W= 0d0
         call eigensystem_c( 'V', 'U', ndim, Ham, W)

         nfill= 0
         do i=1, ndim
            if (W(i).le.Ef)nfill= nfill+ 1
         enddo

         ! > prepare velocity operator
         call velocity(ndim, thetax, thetay, H01x, H01y, vx, vy)
         vx= vx0 + vx
         vy= vy0 + vy

         sigma= zzero
         do i= 1, nfill
            zn= ham(:, i)
            En= w(i)
            do j= nfill+1, Ndim
               zm= ham(:, j)
               Em= w(j)
               call zgemv('N', ndim, ndim, zone, vy, ndim, zm,1, zzero, vec1,1)
               z1= zdotc(ndim, zn, 1, vec1, 1)
               call zgemv('N', ndim, ndim, zone, vx, ndim, zn,1, zzero, vec1,1)
               z2= zdotc(ndim, zm, 1, vec1, 1)
               z3= z1*z2
               call zgemv('N', ndim, ndim, zone, vx, ndim, zm,1, zzero, vec1,1)
               z1= zdotc(ndim, zn, 1, vec1, 1)
               call zgemv('N', ndim, ndim, zone, vy, ndim, zn,1, zzero, vec1,1)
               z2= zdotc(ndim, zm, 1, vec1, 1)
               z3= z3- z1*z2
               z3= z3/(En-Em)**2
               sigma(i)= sigma(i)+ aimag(z3)
            enddo ! j unoccupied
         enddo ! i occupied
         res= res+ sum(sigma(1:nfill))
         print *, ik, res
         write(10, '(100f16.8)')kpoint(ik,:), sum(sigma(1:nfill))
         if (mod(ik,NN).eq.0)write(10, *)''

      enddo ! ik
      res= res/dble(NN*NN)*2d0*pi/dble(Nx*Ny)


      return

   end subroutine chern_kubo

!>>> this method comes from  
! Ref:  Journal of the Physical Society of Japan
! Vol. 73, No. 10, October, 2004, pp. 2624–2627
! Anomalous Hall Effect and Skyrmion Number in Real and Momentum Spaces
! Masaru ONODA1, Gen TATARA2 and Naoto NAGAOSA1;3
! Eq.(3)
! construct on 20th Dec 2012
   subroutine chern_streda(Ef, res)

      use para
      implicit none

      !>> inout variables
      real(dp), intent(in) :: Ef
      real(dp), intent(out) :: res

      integer :: ik, i, j
 
      ! number of occupied states
      integer :: nfill

      integer :: nnzx
      integer :: nnzy
      integer :: nnzmax
 
      real(dp) :: thetax
      real(dp) :: thetay

      real(dp) :: En
      real(dp) :: Em

      real(dp) :: tau

      complex(dp) :: z1
      complex(dp) :: z2
      complex(dp) :: z3

      character :: matdescra(6)

      integer, allocatable :: iax_coo(:)
      integer, allocatable :: jax_coo(:)
      integer, allocatable :: iay_coo(:)
      integer, allocatable :: jay_coo(:)
      complex(dp), allocatable :: vx_coo(:)
      complex(dp), allocatable :: vy_coo(:)

      real(Dp), allocatable :: W(:)
      real(dp), allocatable :: sigma(:)

      ! Hamiltonian of bulk system
      complex(Dp), allocatable :: Ham(:, :) 
      complex(Dp), allocatable :: dis (:,:)
 
      complex(dp), allocatable :: zn(:)
      complex(dp), allocatable :: zm(:)
      complex(dp), allocatable :: vec1(:)

      complex(dp), allocatable ::  X(:, :)
      complex(dp), allocatable ::  Y(:, :)
      complex(dp), allocatable :: vx(:, :)
      complex(dp), allocatable :: vy(:, :)
      complex(dp), allocatable :: vx1(:, :)
      complex(dp), allocatable :: vy1(:, :)
      complex(dp), allocatable :: vx2(:, :)
      complex(dp), allocatable :: vy2(:, :)
      complex(dp), allocatable :: vx0(:, :)
      complex(dp), allocatable :: vy0(:, :)

      complex(dp), allocatable :: H01x(:, :)
      complex(dp), allocatable :: H01y(:, :)

      ! vector multiply vector v1*v2
      complex(dp), external :: zdotc

      allocate(H01x(nband*2, nband*2))
      allocate(H01y(nband*2, nband*2))


      allocate(W(ndim))
      allocate(ham(ndim, ndim))
      allocate(dis(ndim, ndim))
      allocate(zn(ndim))
      allocate(zm(ndim))
      allocate(vec1(ndim))
      allocate(vx(ndim, ndim))
      allocate(vy(ndim, ndim))
      allocate(vx1(ndim, ndim))
      allocate(vy1(ndim, ndim))
      allocate(vx2(ndim, ndim))
      allocate(vy2(ndim, ndim))
      allocate( X(ndim, ndim))
      allocate( Y(ndim, ndim))
      allocate(vx0(ndim, ndim))
      allocate(vy0(ndim, ndim))
      allocate(sigma(ndim))

      nnzmax= 10*ndim
      allocate(iax_coo(nnzmax))
      allocate(iay_coo(nnzmax))
      allocate(jax_coo(nnzmax))
      allocate(jay_coo(nnzmax))
      allocate(vx_coo(nnzmax))
      allocate(vy_coo(nnzmax))
      
      call disorder(dis)

      ! > calculate H01x, H01y
      do ik=1, nrpts
         i = irvec(1, ik)
         j = irvec(2, ik)
         if (i .eq. 1 .and. j .eq. 0)H01x= Hmnr(:, :, ik)
         if (i .eq. 0 .and. j .eq. 1)H01y= Hmnr(:, :, ik)
      enddo

     !call velocity2(ndim, vx0, vy0)

      ! >> prepare eigenvectors
 
      thetax = 0d0
      thetay = 0d0
      tau = 40d0
 


      ! > prepare velocity operator
      Ham=0d0
     !call ham_twist(thetax, thetay, Ham)
     !call ham_bulk(thetax, thetay, Ham)
      call ham_haldane2(thetax, thetay, Ham)
      call ham_xy(X, Y)

      ! vx= -i[X, H]
      call mat_mul(ndim, X, Ham, vx)
      call mat_mul(ndim, Ham, X, vy)
      vx= -zi*(vx- vy)
 
      ! vy= -i[Y, H]
      call mat_mul(ndim, Y, Ham, vy)
      call mat_mul(ndim, Ham, Y, X)
      vy= -zi*(vy- X)

      ! calculation bulk hamiltonian with twist boundary
     !call ham_twist(thetax, thetay, Ham)
      call ham_haldane(thetax, thetay, Ham)
      Ham= Ham+ dis
      W= 0d0
      call eigensystem_c( 'V', 'U', ndim, Ham, W)

      nfill= 0
      do i=1, ndim
         if (W(i).le.Ef)nfill= nfill+ 1
      enddo


      ik= 0
      do i=1, ndim 
         do j=1, ndim
            if (abs(vx(i,j)).ge.1e-3)then 
               ik=ik+ 1
               iax_coo(ik)= i
               jax_coo(ik)= j
                vx_coo(ik)= vx(i,j)
            endif
         enddo
      enddo
      nnzx= ik
      
      ik= 0
      do i=1, ndim 
         do j=1, ndim
            if (abs(vy(i,j)).ge.1e-3)then 
               ik=ik+ 1
               iay_coo(ik)= i
               jay_coo(ik)= j
                vy_coo(ik)= vy(i,j)
            endif
         enddo
      enddo
      nnzy= ik

      matdescra(1)= 'H'
      matdescra(2)= 'U'
      matdescra(3)= 'U'
      matdescra(4)= 'F'

      sigma= zzero
      do i= 1, nfill
         zn= ham(:, i)
         En= w(i)
         do j= nfill+1, Ndim
            zm= ham(:, j)
            Em= w(j)
            !<n|vx|m><m|vy|n>
           !call zgemv('N', ndim, ndim, zone, vy, ndim, zm,1, zzero, vec1,1)
            call mkl_zcoomv('N', ndim, ndim, zone, matdescra, &
                            vy_coo, iay_coo, jay_coo, nnzy, zm, zzero, vec1)
            z1= zdotc(ndim, zn, 1, vec1, 1)
           !call zgemv('N', ndim, ndim, zone, vx, ndim, zn,1, zzero, vec1,1)
            call mkl_zcoomv('N', ndim, ndim, zone, matdescra, &
                            vx_coo, iax_coo, jax_coo, nnzx, zn, zzero, vec1)
            z2= zdotc(ndim, zm, 1, vec1, 1)
            z3= z1*z2

            !<n|vy|m><m|vx|n> is equal to conjugate(<n|vx|m><m|vy|n>)
            !<n|vy|m><m|vx|n>
           !call zgemv('N', ndim, ndim, zone, vx, ndim, zm,1, zzero, vec1,1)
           !call mkl_zcoomv('N', ndim, ndim, zone, matdescra, &
           !                vx_coo, iax_coo, jax_coo, nnzx, zm, zzero, vec1)
           !z1= zdotc(ndim, zn, 1, vec1, 1)
           !call zgemv('N', ndim, ndim, zone, vy, ndim, zn,1, zzero, vec1,1)
           !call mkl_zcoomv('N', ndim, ndim, zone, matdescra, &
           !                vy_coo, iay_coo, jay_coo, nnzy, zn, zzero, vec1)
           !z2= zdotc(ndim, zm, 1, vec1, 1)

            !<n|vx|m><m|vy|n>-<n|vy|m><m|vx|n>
           !z3= z3- z1*z2
            z3= z3*tau*tau/(1d0+(En-Em)**2*tau*tau)
            sigma(i)= sigma(i)+ aimag(z3)
         enddo ! j unoccupied
      enddo ! i occupied
      res= res+ sum(sigma(1:nfill))*2

      res= res/dble(Nx*Ny)*2d0*pi


      return

   end subroutine chern_streda


   ! this subroutine is wrong
   ! if H= Cn^+ C1+ C1^+ Cn
   ! then vx= -i [X, H]= (1-n)C1^+ Cn -(1-n)Cn^+C1
   ! but we need not this subroutine any more
   ! for comparasion, we keep this
   ! Dec 25th 2012, the first Christmas day after the end of the world
   subroutine velocity(ndim, thetax, thetay, H01x, H01y, vx, vy)

      use para,only : dp, Nband, zzero, zi, Nx, Ny

      implicit none

      !>> inout variables
      integer, intent(in) :: ndim

      real(dp), intent(in) :: thetax
      real(dp), intent(in) :: thetay

      complex(dp), intent(in) :: H01x(Nband*2, Nband*2)
      complex(dp), intent(in) :: H01y(Nband*2, Nband*2)

      complex(dp), intent(out) :: vx(ndim, ndim)
      complex(dp), intent(out) :: vy(ndim, ndim)

      !>> local variables

      integer :: ix1
      integer :: ix2
      integer :: iy
      integer :: idx1
      integer :: idx2

      complex(dp) :: ratio

      vx= zzero
      vy= zzero

      !> vx
      ix1= Nx
      ix2= 1

      ratio= cos(thetax)+ zi*sin(thetax)
      do iy=1, Ny
         idx1= (ix1-1)*Ny*Nband*2+(iy-1)*Nband*2
         idx2= (ix2-1)*Ny*Nband*2+(iy-1)*Nband*2
         vx(idx1+1: idx1+Nband*2, idx2+1: idx2+Nband*2)= &
         vx(idx1+1: idx1+Nband*2, idx2+1: idx2+Nband*2)+ &
         zi*ratio*H01x
      enddo

      ix1= 1
      ix2= Nx

      do iy=1, Ny
         idx1= (ix1-1)*Ny*Nband*2+(iy-1)*Nband*2
         idx2= (ix2-1)*Ny*Nband*2+(iy-1)*Nband*2
         vx(idx1+1: idx1+Nband*2, idx2+1: idx2+Nband*2)= &
         vx(idx1+1: idx1+Nband*2, idx2+1: idx2+Nband*2)- &
         zi*conjg(ratio)*conjg(transpose(H01x))
      enddo

      !> vy
      ix1= Ny
      ix2= 1

      ratio= cos(thetay)+ zi*sin(thetay)
      do iy=1, Nx
         idx1= (iy-1)*Ny*Nband*2+(ix1-1)*Nband*2
         idx2= (iy-1)*Ny*Nband*2+(ix2-1)*Nband*2
         vy(idx1+1: idx1+Nband*2, idx2+1: idx2+Nband*2)= &
         vy(idx1+1: idx1+Nband*2, idx2+1: idx2+Nband*2)+ &
         zi*ratio*H01y
      enddo

      ix1= 1
      ix2= Ny

      do iy=1, Nx
         idx1= (iy-1)*Ny*Nband*2+(ix1-1)*Nband*2
         idx2= (iy-1)*Ny*Nband*2+(ix2-1)*Nband*2
         vy(idx1+1: idx1+Nband*2, idx2+1: idx2+Nband*2)= &
         vy(idx1+1: idx1+Nband*2, idx2+1: idx2+Nband*2)- &
         zi*conjg(ratio)*conjg(transpose(H01y))
      enddo

      return

   end subroutine velocity

   subroutine velocity2(ndim, vx, vy)

      use para,only : dp, Nband, zzero, zi, Nx, Ny, ijmax

      implicit none

      !>> inout variables
      integer, intent(in) :: ndim

      complex(dp), intent(out) :: vx(ndim, ndim)
      complex(dp), intent(out) :: vy(ndim, ndim)

      !>> local variables

      ! loop index  
      integer :: i1
      integer :: i2
      integer :: j1
      integer :: j2
      
      complex(Dp), allocatable :: Hij(:, :, :, :)
      
      allocate(Hij(-ijmax:ijmax, -ijmax:ijmax, nband*2, nband*2))
      
      call ham_qlayer2qlayer3(Hij)
      
      
      vx=0.0d0 
      do i1= 1, Nx
      do j1= 1, Ny
      do i2= 1, Nx   
        if (i2.gt.i1.and.(i2-i1).le.ijmax) then
            vx((i1-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                       (i1-1)*Nband*2*Ny+j1*Nband*2,       &
                       (i2-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                       (i2-1)*Nband*2*Ny+j1*Nband*2 )      &
            = zi*Hij(i2-i1, 0, 1:nband*2, 1:nband*2)
        endif 
        if (i1.gt.i2.and.(i1-i2).le.ijmax) then
            vx((i1-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                       (i1-1)*Nband*2*Ny+j1*Nband*2,       &
                       (i2-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                       (i2-1)*Nband*2*Ny+j1*Nband*2 )      &
            = -zi*Hij(i2-i1, 0, 1:nband*2, 1:nband*2)
        endif 
      enddo
      enddo
      enddo

      
      vy=0.0d0 
      do i1= 1, Nx
      do j1= 1, Ny
      do j2= 1, Ny   
        if (j2.gt.j1.and.(j2-j1).le.ijmax) then
            vy((i1-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                       (i1-1)*Nband*2*Ny+j1*Nband*2,       &
                       (i1-1)*Nband*2*Ny+(j2-1)*Nband*2+1: &
                       (i1-1)*Nband*2*Ny+j2*Nband*2 )      &
            = zi*Hij(0, j2-j1, 1:nband*2, 1:nband*2)
        endif 
        if (j1.gt.j2.and.(j1-j2).le.ijmax) then
            vy((i1-1)*Nband*2*Ny+(j1-1)*Nband*2+1: &
                       (i1-1)*Nband*2*Ny+j1*Nband*2,       &
                       (i1-1)*Nband*2*Ny+(j2-1)*Nband*2+1: &
                       (i1-1)*Nband*2*Ny+j2*Nband*2 )      &
            = -zi*Hij(0, j2-j1, 1:nband*2, 1:nband*2)
        endif 
      enddo
      enddo
      enddo

! check hermitcity

    !do i1=1, Ndim
    !do i2=i1+1, Ndim
    !   if(abs(vx(i1,i2)-conjg(vx(i2,i1))).ge.1e-6)then
    !     write(*,*)'there are something wrong with vx'
    !     stop
    !   endif 
    !   if(abs(vy(i1,i2)-conjg(vy(i2,i1))).ge.1e-6)then
    !     write(*,*)'there are something wrong with vy'
    !     stop
    !   endif 
    !enddo
    !enddo

      return

   end subroutine velocity2
