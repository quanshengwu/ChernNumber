!+---------+---------+---------+---------+---------+---------+--------+!
! This subroutine is used for providing disorder H00 and H01 
! Disorder type : 1. Simple  2. magnetic  3. Vacuum  4. only on surface
! Simple type     : only change onsite energy randomly with spin equally
! Magnetic type   : change onsite energy randomly with spin unequally
! Vacuum type     : delete some atoms randomly
! Constructed by Quan Sheng Wu on Jau/8/2011
!+---------+---------+---------+---------+---------+---------+--------+!

   subroutine disorder(dh00) 

      use para
      implicit none

      ! loop index
      integer :: i
      integer :: i1
      integer :: n1

      real(Dp) :: dh, dh1

      ! Hamiltonian in principle layer
      complex(Dp),intent(inout) :: dh00(Ndim,Ndim)

      dh00=zzero

      IF (disorder_type .eq. 'simple')THEN

      ! Simple type
      do n1=1, Ny ! y direction
      do i1=1, Nx ! x direction

         do i= 1, Nband
            ! disorder rate
            if(ran(Seed).le.rate)then
            ! spin up
            dh=(-0.5d0+ran(Seed))*disorder_strength
            dh00((i1-1)*Ny*Nband*2+(n1-1)*Nband*2+i,&
            (i1-1)*Ny*Nband*2+(n1-1)*Nband*2+i)=dh
            ! spin down 
            dh=(-0.5d0+ran(Seed))*disorder_strength
            dh00((i1-1)*Ny*Nband*2+(n1-1)*Nband*2+i+Nband,&
            (i1-1)*Ny*Nband*2+(n1-1)*Nband*2+i+Nband)=dh
            endif
         enddo
 
      enddo!i1
      enddo!n1
      
      ELSEIF (disorder_type .eq. 'surface')THEN

      ! only on surface
      do n1=1,Ny ! y direction
      do i1=1,Nx ! x direction

         ! diagnoal
         ! on the surface
         dh=(-0.5d0+ran(Seed))*disorder_strength
         if (i1.eq.1 .or. i1 .eq. Ny)then 
            do i=1,Nband
            ! disorder rate
            if (ran(Seed).le.rate)then
               !          print*, 'there are some random impurites'
               ! spin up
               dh00((i1-1)*Ny*Nband*2+(n1-1)*Nband*2+i,&
               (i1-1)*Ny*Nband*2+(n1-1)*Nband*2+i)=dh
               ! spin down 
               dh00((i1-1)*Ny*Nband*2+(n1-1)*Nband*2+i+Nband,&
               (i1-1)*Ny*Nband*2+(n1-1)*Nband*2+i+Nband)=dh
            endif
            enddo
         endif
  
      enddo!i1
      enddo!n1

      ELSEIF (disorder_type .eq. 'bulk_magnetic')THEN

      ! bulk magnetic impurity type
      do i1=1,Ny ! y direction
      do n1=1,Nx ! x direction

         do i=1, Nband
         ! disorder rate
         if(ran(Seed).le.rate)then
         !          print*, 'there are some random impurites'
         dh =(-0.5d0+ran(Seed))*disorder_strength
         dh1=(-0.5d0+ran(Seed))*disorder_strength
         ! spin up
         dh00((i1-1)*Np*Nband*2+(n1-1)*Nband*2+i,&
         (i1-1)*Np*Nband*2+(n1-1)*Nband*2+i)=dh
         ! spin down 
         dh00((i1-1)*Np*Nband*2+(n1-1)*Nband*2+i+Nband,&
         (i1-1)*Np*Nband*2+(n1-1)*Nband*2+i+Nband)=dh1
         endif
         enddo
 
      enddo!i1
      enddo!n1
      
      ELSE

          Write(*,*)' no disorder type input'

      ENDIF


      return

   end subroutine disorder
