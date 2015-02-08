! > haldane's model
   subroutine readHmnR()

      use para
      implicit none

      integer :: ir

      irvec=0
      ndegen=1
      hmnr= zzero

      ! 0 0 0
      ir= 1
      irvec(1, ir)= 0
      irvec(2, ir)= 0
      hmnr(1, 1, ir)= M- 4d0*B
      hmnr(2, 2, ir)=-M+ 4d0*B

      ! 1 0 
      ir= 2
      irvec(1, ir)= 1
      irvec(2, ir)= 0
      hmnr(1, 1, ir)= B
      hmnr(2, 2, ir)=-B
      hmnr(1, 2, ir)= zi*A/2d0
      hmnr(2, 1, ir)= zi*A/2d0 

      ! 0 1 
      ir= 3
      irvec(1, ir)= 0
      irvec(2, ir)= 1
      hmnr(1, 1, ir)= B
      hmnr(2, 2, ir)=-B
      hmnr(1, 2, ir)=-A/2d0
      hmnr(2, 1, ir)= A/2d0


      !-1 0 
      ir= 4
      irvec(1, ir)=-1
      irvec(2, ir)= 0
      hmnr(1, 1, ir)= B
      hmnr(2, 2, ir)= -B
      hmnr(1, 2, ir)= -zi*A/2d0
      hmnr(2, 1, ir)= -zi*A/2d0 

      ! 0-1 
      ir= 5
      irvec(1, ir)= 0
      irvec(2, ir)=-1
      hmnr(1, 1, ir)= B
      hmnr(2, 2, ir)=-B
      hmnr(1, 2, ir)= A/2d0
      hmnr(2, 1, ir)=-A/2d0



      return
   end subroutine
