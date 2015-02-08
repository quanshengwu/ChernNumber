
   subroutine chern

      use mpi
      use para
      implicit none

      integer :: i
      integer :: j
      integer :: idis
      integer :: ierr

      real(dp) :: Ef

      real(dp) :: res

      real(dp), allocatable :: omega(:)
      real(dp), allocatable :: chernnumber(:, :)
      real(dp), allocatable :: chernnumber_mpi(:, :)

      allocate(omega(omeganum))
      allocate(chernnumber(omeganum, Ndis))
      allocate(chernnumber_mpi(omeganum, Ndis))
      chernnumber_mpi= 0d0

      do i=1, omeganum
         if (omeganum.ne.1)then
            omega(i)=minomega+(i-1)*(maxomega-minomega)/(omeganum-1)
         else
            omega(i)=minomega
         endif
      enddo

      ! > using mpi
      do idis= 1+cpuid, Ndis, num_cpu
         Seed= Seed+ idis
         do i=1, omeganum
            Ef= omega(i)
           !call chern_realspace( Ef, res)
           !call chern_prodan( Ef, res)
           !call chern_kubo( Ef, res)
            call chern_streda( Ef, res)
            chernnumber(i, idis)= res
            if (cpuid.eq.0)write(*, '(i,2f16.8)')idis, Ef, res
         enddo
      enddo

      call mpi_allreduce(chernnumber, chernnumber_mpi, size(chernnumber), &
                         mpi_dp, mpi_sum, mpi_cw, ierr)
      chernnumber= chernnumber_mpi

      if (cpuid.eq.0) then
         open(unit=20, file=filename)
         do i=1, Omeganum
            write(20, '(100000f16.8)')omega(i), sum(chernnumber(i, :))/dble(Ndis), &
                                    chernnumber(i, :)
         enddo
         close(20)
      endif

      return
   end subroutine chern
