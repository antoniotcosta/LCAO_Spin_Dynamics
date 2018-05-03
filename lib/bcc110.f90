subroutine bcc110

 use f90_kind
 use constants_u
 use lattice

 real(double) :: aux

!         n01,n02   = # of 1st. and 2nd. in-plane n.n. respectively
!         n1,n2     = # of 1st. and 2nd. inter-plane n.n. respectively
!         r0( i, l) = coord. of the in-plane (first & second n.n.)
!                     (i=1,n0; l = 1,3 (x,y,z) in units of (a0/2);
!                     a0 = latt. const. 
!         c0( i, l) = direction cosines of r0
!         r1( i, l) = coord. of the inter-plane 1st. n.n. (i=1,n1)
!         c1( i, l) = direction cosines of r1
!         r2( i, l) = coord. of the inter-plane 2nd. n.n. (i=1,n2)
!         c2( i, l) = direction cosines of r2 

! in-plane 1st. and 2nd. n.n.

  n01 = 4
  n02 = 2
      
  r0(1,1) = 1.d0
  r0(1,2) = 1.d0
  r0(1,3) = 1.d0

  r0(2,1) =-1.d0
  r0(2,2) =-1.d0
  r0(2,3) = 1.d0

  r0(3,1) =-1.d0
  r0(3,2) =-1.d0
  r0(3,3) =-1.d0

  r0(4,1) = 1.d0
  r0(4,2) = 1.d0
  r0(4,3) =-1.d0

  r0(5,1) = 0.d0
  r0(5,2) = 0.d0
  r0(5,3) = 2.d0

  r0(6,1) = 0.d0
  r0(6,2) = 0.d0
  r0(6,3) =-2.d0

  aux = 1.d0/sq3
      
  c0(1,1) = aux
  c0(1,2) = aux
  c0(1,3) = aux

  c0(2,1) =-aux
  c0(2,2) =-aux
  c0(2,3) = aux

  c0(3,1) =-aux
  c0(3,2) =-aux
  c0(3,3) =-aux

  c0(4,1) = aux
  c0(4,2) = aux
  c0(4,3) =-aux 

  c0(5,1) = 0.d0
  c0(5,2) = 0.d0
  c0(5,3) = 1.d0

  c0(6,1) = 0.d0
  c0(6,2) = 0.d0
  c0(6,3) =-1.d0

!  inter-plane 1st. n.n.

  n1=2

  r1(1,1) =-1.0d0
  r1(1,2) = 1.0d0
  r1(1,3) = 1.0d0

  r1(2,1) =-1.0d0
  r1(2,2) = 1.0d0
  r1(2,3) =-1.0d0

  c1(1,1) =-aux
  c1(1,2) = aux
  c1(1,3) = aux

  c1(2,1) =-aux
  c1(2,2) = aux
  c1(2,3) =-aux

!  inter-plane 2nd. n.n.

  n2=2

  r2(1,1) = 0.d0
  r2(1,2) = 2.d0
  r2(1,3) = 0.d0

  r2(2,1) =-2.d0
  r2(2,2) = 0.d0
  r2(2,3) = 0.d0

  c2(1,1) = 0.d0
  c2(1,2) = 1.d0
  c2(1,3) = 0.d0

  c2(2,1) =-1.d0
  c2(2,2) = 0.d0
  c2(2,3) = 0.d0
       
  return

end subroutine bcc110      
