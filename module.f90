  module para
!    use mpi

     implicit none
     include 'mpif.h'

! output file name
     character*40 :: filename
    
! disorder type 
     character*40 :: disorder_type

! some constant for mpi 
     integer,parameter :: mpi_dp=mpi_double_precision
     integer,parameter :: mpi_dc=mpi_double_complex
     integer,parameter :: mpi_cw=mpi_comm_world
     integer,parameter :: mpi_in=mpi_integer

! double precision  
     integer,parameter :: Dp=kind(1.0d0)

! number of bands in a primitive cell, not include spin degenercy 
     integer,parameter :: Nband=1

     integer, parameter :: ijmax=2

! Nx is the length of nanoribbon along[11-20] direction,
! we set x direction the growing direction, Nx is the width 
! of the centre zone.
     integer :: Nx

! Ny is the width of nanoribbon 
     integer :: Ny

! Np is principle layer thickness along x direction
     integer :: Np

! Ndim is the dimension of H00, usually Ndim=30*Ny*Np 
! 30 is the total orbitals in a unit cell
     integer :: ndim

! number of R points, read from HmnR.data
     integer :: Nrpts

! a parameter to control soc
! Soc=0 means no spin-orbit coupling
! Soc!=0 means no spin-orbit coupling
     integer :: Soc

! cpu index
     integer :: cpuid

! cpu number we used 
     integer :: num_cpu

! number of gpu
     integer :: num_gpu

! the number of B 
     integer :: numB 

! the number of disorder samples
     integer :: Ndis

! the number of omega
     integer :: omeganum 

!  random number seed
     integer :: Seed

! Gpu device number
     integer :: gpu_device

! disorder rate
     real(Dp) :: Rate

! disorder strength , onsite energy maximum amplitude
     real(Dp) :: disorder_strength

! a infinite small number 
     real(Dp) :: eta 

! magnetic strength interval 
     real(Dp) ::  minB
     real(Dp) ::  maxB

! omega interval 
     real(Dp) ::  minomega
     real(Dp) ::  maxomega

! circumference ratio pi  
     real(dp),parameter :: Pi= 3.14159265359d0

! complex constant 
     complex(dp),parameter    :: zi=(0.0d0, 1.0d0)
     complex(dp),parameter    :: zone=(1.0d0, 0.0d0)
     complex(dp),parameter    :: zzero=(0.0d0, 0.0d0)

! parameters for H(k)

     ! lattice constant
     real(dp), parameter :: M  =  0.40d0
     real(dp), parameter :: phi=  Pi/3d0
     real(dp), parameter :: t1 =  0.20d0
     real(dp), parameter :: t0 =  1.00d0

     real(dp), parameter :: A = 1d0
     real(dp), parameter :: AA= 1d0
     real(dp), parameter :: B =-1d0
     real(dp), parameter :: C = 0d0
     real(dp), parameter :: D = 0d0
    !real(dp), parameter :: M =-2d0

     complex(dp), parameter :: vsp= -zi*A/2d0/AA
     complex(dp), parameter :: vss= (B+D)/AA/AA
     complex(dp), parameter :: vpp= (-B+D)/AA/AA
     complex(dp), parameter :: Es= C+M-4d0*(B+D)/AA/AA
     complex(dp), parameter :: Ep= C-M-4d0*(-B+D)/AA/AA
        
        
! R coordinates  
     integer, allocatable     :: Irvec(:,:)

! Hamiltonian m,n are band indexes
     complex(dp), allocatable :: HmnR(:,:,:)

! no of degeneracy of R point 
     integer, allocatable     :: ndegen(:)

     contains

     subroutine now(time_now)
     
        implicit none
     
        integer   :: time_new(8)
        real(Dp)  :: time_now
     
     
        call Date_and_time(values=time_new)
     
        time_now=time_new(3)*24*3600+time_new(5)*3600+&
                 time_new(6)*60+time_new(7)+time_new(8)/1000d0
     
        return
     end subroutine now
     
     subroutine show_now
     
        implicit none
        integer   :: time_new(8)
     
        call Date_and_time(values=time_new)
     
        write(*,"(a11,i4,a,i2,a,i2,a2,i2,a,i2,a,i2)")&
        '   Now is: ',time_new(1),'/',time_new(2),'/',time_new(3),' ',&
        time_new(5),':',time_new(6),':',time_new(7)
     
        return
     end subroutine show_now

 end module para
