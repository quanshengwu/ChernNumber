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
      hmnr(1, 1, ir)=  m 
      hmnr(2, 2, ir)= -m
      hmnr(1, 2, ir)=  t0
      hmnr(2, 1, ir)=  t0

      ! 1 0 
      ir= 2
      irvec(1, ir)= 1
      irvec(2, ir)= 0
      hmnr(1, 1, ir)= (cos(phi)-zi*sin(phi)) *t1
      hmnr(2, 2, ir)= (cos(phi)+zi*sin(phi)) *t1
      hmnr(1, 2, ir)= t0

      ! 0 1 
      ir= 3
      irvec(1, ir)= 0
      irvec(2, ir)= 1
      hmnr(1, 1, ir)= (cos(phi)-zi*sin(phi)) *t1
      hmnr(2, 2, ir)= (cos(phi)+zi*sin(phi)) *t1
      hmnr(2, 1, ir)= t0

      ! 1 1 
      ir= 4
      irvec(1, ir)= 1
      irvec(2, ir)= 1
      hmnr(1, 1, ir)= (cos(phi)+zi*sin(phi)) *t1
      hmnr(2, 2, ir)= (cos(phi)-zi*sin(phi)) *t1

      !-1 0 
      ir= 5
      irvec(1, ir)=-1
      irvec(2, ir)= 0
      hmnr(1, 1, ir)= (cos(phi)+zi*sin(phi)) *t1
      hmnr(2, 2, ir)= (cos(phi)-zi*sin(phi)) *t1
      hmnr(2, 1, ir)= t0

      ! 0-1 
      ir= 6
      irvec(1, ir)= 0
      irvec(2, ir)=-1
      hmnr(1, 1, ir)= (cos(phi)+zi*sin(phi)) *t1
      hmnr(2, 2, ir)= (cos(phi)-zi*sin(phi)) *t1
      hmnr(1, 2, ir)= t0

      !-1-1 
      ir= 7
      irvec(1, ir)=-1
      irvec(2, ir)=-1
      hmnr(1, 1, ir)= (cos(phi)-zi*sin(phi)) *t1
      hmnr(2, 2, ir)= (cos(phi)+zi*sin(phi)) *t1

      return
   end subroutine
