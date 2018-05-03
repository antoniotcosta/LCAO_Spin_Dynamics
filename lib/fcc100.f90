subroutine fcc100 

 use f90_kind
 use constants_u
 use lattice

 integer :: i

!   in plane:
  n01 = 4
  n02 = 4

  r0(1,1) = 1.D0
  r0(1,2) = 1.D0
  r0(1,3) = 0.d0

  r0(2,1) = 1.D0
  r0(2,2) = -1.D0
  r0(2,3) =  0.d0

  r0(3,1) = -1.D0
  r0(3,2) = 1.D0
  r0(3,3) =  0.d0

  r0(4,1) = -1.D0
  r0(4,2) = -1.D0
  r0(4,3) =  0.d0

  r0(5,1) = 2.d0
  r0(5,2) = 0.d0
  r0(5,3) = 0.d0

  r0(6,1) = 0.d0
  r0(6,2) = 2.d0
  r0(6,3) = 0.d0
  
  r0(7,1) = 0.d0
  r0(7,2) = -2.d0
  r0(7,3) = 0.d0

  r0(8,1) = -2.d0
  r0(8,2) = 0.d0
  r0(8,3) = 0.d0

  do i=1,n01
   c0(i,:) = r0(i,:)/(sq2)
  end do

  do i=n01+1,n01+n02
   c0(i,:) = r0(i,:)/(2.d0)
  end do

!  inter plane
  n1 = 4
  n2 = 1

  r1(1,1) = 1.D0
  r1(1,2) = 0.d0
  r1(1,3) = 1.D0

  r1(2,1) = 0.d0
  r1(2,2) = 1.D0
  r1(2,3) = 1.D0

  r1(3,1) =  0.d0
  r1(3,2) = -1.D0
  r1(3,3) =  1.D0

  r1(4,1) = -1.D0
  r1(4,2) =  0.d0
  r1(4,3) =  1.D0

  do i=1,4
   r1(i+4,1:2) = r1(i,1:2)
   r1(i+4,3) = -r1(i,3)
  end do

  r2(1,1) = 0.d0
  r2(1,2) = 0.d0
  r2(1,3) = 2.d0

  r2(2,1) = 0.d0
  r2(2,2) = 0.d0
  r2(2,3) = -2.d0

  c1 = r1/sq2
  c2 = r2/(2.d0)
       
  return
 
end
