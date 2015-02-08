! this subroutine is used to read some paramters from
! input.dat
! constructed on Dec/25/2010 by QS.Wu 


  subroutine readinput

     use para
     implicit none

     character*80 :: fname='input.dat'
     logical ::  exists
    
     INQUIRE (FILE = fname, EXIST = exists)
 
     inquire(file=fname,exist=exists)
     if (exists)then
        write(*,*) '  > Read some paramters from input.dat <'
        open(unit=1001,file=fname,status='old')
     else
        write(*,*)'File ' ,fname, 'dosnot exist'
        stop
     endif

     read(1001,*)filename
     write(*,'(2x,a16,11x,a,a42)')'Output file name ',':',filename
     read(1001,*)gpu_device
     write(*,'(2x,a14,11x,a,i)')'Gpu device id ',':',gpu_device
     read(1001,*)omeganum,numB
     write(*,'(2x,a13,14x,a,i5,a20,i5)')'Omega number ',&
                                        ':',omeganum,'B number: ',numB
     read(1001,*)minomega,maxomega
     write(*,'(2x,a15,12x,a,2f9.5)')'Omega interval ',':',minomega,maxomega
     read(1001,*)minB,maxB
     write(*,'(2x,a27,a,2f9.3)')'Magnetic strength interval ',':',minB,maxB
     read(1001,*)disorder_type
     write(*,'(2x,a14,13x,a,a20)')'Disorder type ',':',disorder_type
     read(1001,*)Rate,disorder_strength
     write(*,'(2x,a14,13x,a,f10.3)')'Disorder rate ',':',Rate
     write(*,'(2x,a20,7x,a,f10.3)')'Disorder strength ',':',disorder_strength
     read(1001,*)Ndis
     write(*,'(2x,a17,10x,a,i10)')'Number of disorder ',':',Ndis
     read(1001,*)Seed
     write(*,'(2x,a17,10x,a,i10)')'Random number seed',':',Seed
     read(1001,*)Nx
     write(*,'(2x,a3,24x,a,i10)')'Nx ',':',Nx
     read(1001,*)Ny,Np
     write(*,'(2x,a9,18x,a,3i5)')'Ny,Np ',':',Ny,Np
     read(1001,*)Soc
     write(*,'(2x,a9,18x,a,i5)')'Soc flag ',':',Soc
     read(1001,*)eta
     write(*,'(2x,a4,23x,a,E18.4)')'Eta ',':',eta

     Ndim=Nband*2*Ny*Nx
 
     write(*,'(2x,a4 ,23x,a,i5)')'Ndim',':',Ndim


     close(1001)

     write(*,*)'  >> Read input.dat file successfully <<'
     return
  end subroutine



