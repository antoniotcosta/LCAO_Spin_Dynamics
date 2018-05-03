subroutine bcc100

 use f90_kind
 use constants_u
 use lattice

 real(double) :: aux0

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

!     in common:
!         r0( i, l) = coord.  of  the  in-plane nearest-neighbours (i=1,n0)
!                     l = 1,3 (x,y,z) in units of (a0/2); a0 = latt. const.
!         r1( i, l) = coord. of the inter-plane nearest-neighbours (i=1.n1)
!                     l = 1,3 (x,y,z) in units of (a0/2); a0 = latt. const.
!         r2( i, l) = coord. of the in-plane 2nd. n.n. (i=1,n2), (l=1,3) in
!                     units of a0/2 (latt. const.)
!         r3        = coord. of the inter-plane 2nd. n.n. (i=1,n3), idem,
!                     idem
!         cn( i, l) = direction cosines of rn

! in-plane 1st. and 2nd. n.n.

 n01=0
 n02=4

 r0(1,:) = (/ 2.d0, 0.d0, 0.d0 /)

 r0(2,:) = (/ 0.d0, 2.d0, 0.d0 /)

 r0(3,:) = (/ -2.d0, 0.d0, 0.d0 /)

 r0(4,:) = (/ 0.d0, -2.d0, 0.d0 /)

 c0 = r0/2.d0


!  inter-plane 1st. n.n.

 n1=4

 r1(1,:) = (/ 1.0d0, 1.d0, 1.d0 /)

 r1(2,:) = (/ -1.0d0, 1.d0, 1.d0 /)

 r1(3,:) = (/ -1.0d0, -1.0d0, 1.0d0 /)

 r1(4,:) = (/ 1.0d0, -1.0d0, 1.0d0 /)

 aux0 = 1.0d0/sq3

 c1 = r1*aux0

! inter-plane 2nd. n.n.

 n2=1

 r2(1,:) = (/ 0.d0, 0.d0, 2.d0/)

 c2 = r2/2.d0

 return

end subroutine bcc100
