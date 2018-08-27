!--------+--------+--------+--------+--------+--------+--------+------!
! main program for electronic transport properties of 3D topological 
! insulator nanoribbon, using non-equibrium green's function
! constructed by Q.S.Wu on Dec/25/2010
! mail: wuquansheng@gmail.com
!       wuq@phys.ethz.ch
!       spring10boy@163.com
!--------+--------+--------+--------+--------+--------+--------+------!

  program main
    !use mpi
     use para

     implicit none
     
     character *20     :: cpu_nam
     integer           :: namelen 
     

     ! an error index, if ierr !=0 then error occurs
     integer  :: ierr
     integer  :: i
     
     ! Fermi energy
     real(Dp),allocatable :: omega(:)

     ! conductance vs Fermi energy 
     real(Dp),allocatable :: T(:, :)
     real(Dp),allocatable :: T_mpi(:, :)

     ! time measure
     real(Dp) :: time1,time2,time3

     ! initial the environment of mpi 
     call MPI_INIT(ierr)
     call MPI_Comm_rank(MPI_comm_world,cpuid,ierr)
     call mpi_comm_size(MPI_comm_world,num_cpu,ierr)
     call mpi_get_processor_name(cpu_nam,namelen,ierr)


     if (cpuid.eq.0)then
        write(*,*)'  +-----+-----+-----+-----+-----+-----+-----+-----+----+'
        write(*,*)'  +            Begin our program JUPITER               +'
        write(*,*)'  +-----+-----+-----+-----+-----+-----+-----+-----+----+'
        call show_now
     endif

     ! if mpi initial wrong, alarm 
     if(cpuid.eq.0)then
     if(ierr.ne.0)then
       write(*,*)'  >>> Error : mpi initialize wrong'
       stop 
     endif
     endif

     call now(time2)
     time1=time2

     !read some parameters from input.dat
     if(cpuid.eq.0)then
       call readinput
     endif

     nrpts = 7
     allocate(ndegen(nrpts))
     allocate(irvec(2,nrpts))
     allocate(HmnR(nband*2,nband*2,nrpts))

     if (cpuid.eq.0) write(*,*) '  > Begin reading Hmn_R.data <'
     call readHmnR() 
     if (cpuid.eq.0) write(*,*) '  >> Read Hmn_R.data successfully <<'

     ! broadcast nband and Nrpts to every cpu
     call MPI_bcast(Nrpts,1,mpi_in,0,mpi_cw,ierr)

     ! broadcast ndim,Nk,omeganum,Maxomega,nslab,soc,eta to every cpu
     call MPI_bcast(Ndim,1,mpi_in,0,mpi_cw,ierr)
     call MPI_bcast(Nx,1,mpi_in,0,mpi_cw,ierr)
     call MPI_bcast(Ny,1,mpi_in,0,mpi_cw,ierr)
     call MPI_bcast(Np,1,mpi_in,0,mpi_cw,ierr)
     call MPI_bcast(Seed,1,mpi_in,0,mpi_cw,ierr)
     call MPI_bcast(omeganum,1,mpi_in,0,mpi_cw,ierr)
     call MPI_bcast(numB,1,mpi_in,0,mpi_cw,ierr)
     call MPI_bcast(ndis,1,mpi_in,0,mpi_cw,ierr)
     call MPI_bcast(soc,1,mpi_in,0,mpi_cw,ierr)
     call MPI_bcast(eta,1,mpi_dp,0,mpi_cw,ierr)
     call MPI_bcast(minB,1,mpi_dp,0,mpi_cw,ierr)
     call MPI_bcast(maxB,1,mpi_dp,0,mpi_cw,ierr)
     call MPI_bcast(minomega,1,mpi_dp,0,mpi_cw,ierr)
     call MPI_bcast(maxomega,1,mpi_dp,0,mpi_cw,ierr)
     call MPI_bcast(Rate,1,mpi_dp,0,mpi_cw,ierr)
     call MPI_bcast(disorder_type,40,mpi_char,0,mpi_cw,ierr)
     call MPI_bcast(disorder_strength,1,mpi_dp,0,mpi_cw,ierr)

     allocate(omega(omeganum))
     allocate(T(omeganum,Ndis))
     allocate(T_mpi(omeganum,Ndis))
     omega= 0d0
     T = 0d0
     T_mpi= 0d0

     do i=1,omeganum
        if (omeganum.ne.1)then
           omega(i)=minomega+(i-1)*(maxomega-minomega)/(omeganum-1)
        else
           omega(i)=minomega
        endif
     enddo

    !call ek_bulk
     call chern


101  format(3i5,a10,f9.4,a11,f5.2,a3,f6.1,a2,2f10.5)
 
     if(cpuid.eq.0)write( *,*)'  >> End calculate conductance <<'

     ! gather conductance T from each cpu and write it to a file
    !call mpi_allreduce(T, T_mpi, size(T), mpi_dp, mpi_sum, mpi_cw, ierr)
    !T= T_mpi


     if (cpuid.eq.0)then
        do i=1,omeganum
        enddo
     endif

     call now(time3)
     if (cpuid.eq.0)then
        write(*,*)' '
        write(*,*)' +-----+-----+-----+-----+-----+-----+-----+-----+----+'
        write(*,'(a42,f27.1,a2,a5)')'  +>>>  Costing time for all this program :&
        ',time3-time1,' s',' <<<+' 
        call show_now
        write(*,*)' +>> Congratulations! you finished the calculation. <<+'
        write(*,*)' +-----+-----+-----+-----+-----+-----+-----+-----+----+'
     endif
     call mpi_finalize(ierr)

# if defined (__CUDA__)
     ! shutdown cula framework
     call cula_shutdown()
     call cuda_destroy() ! destroy cuda
     call cuda_deallocate
# endif

  end program
