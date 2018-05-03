subroutine helphfcc(ax,ay,h00r,h00i,h01r,h01i,h02r,h02i,mm)

 use f90_kind
 use constants_u
 use lattice
 use tight_binding

 integer :: i,inn,mm
 real(double), dimension(9,9) :: h00r,h00i,h01r,h01i,h02r,h02i,bp
 real(double) :: w(3),ax,ay,spkr0,ckr0,skr0,spkr1,ckr1,skr1,spkr2,ckr2,skr2, &
                 g1,g2,g3,g4,g5,g6,g7,g8,g9,g10


!     input:
!         ax, ay = components of the wave vector k in units of 
!                      2*pi/a0 (a0 = lattice constant)
!         m = 1 we obtain the hamltonian for W
!           = 2 we obtain the hamiltonian for paramagnetic Fe
!
!     in common:
!         pi, sq2 and sq3,
!         n01,n02   = # of 1st. & 2nd. in-plane n.n. respectively
!         n1,n2     = # of 1st. & 2nd. inter-plane n.n. respectively
!         r0( i, l) = coord. of the in-plane (first & second n.n.)
!                     (i=1,n0; l = 1,3 (x,y,z) in units of (a0/2);
!                     a0 = latt. const. 
!         c0( i, l) = direction cosines of r0
!         r1( i, l) = coord. of the inter-plane 1st. n.n. (i=1,n1)
!         c1( i, l) = direction cosines of r1
!         r2( i, l) = coord. of the inter-plane 2nd. n.n. (i=1,n2)
!         c2( i, l) = direction cosines of r2 
! 
!     output:
!         h00r = real part of q-parallel in-plane hoppings 
!         h00i = im.  part of q-parallel in-plane hoppings 
!         h01r = real part of q-parallel inter-plane hoppings (1st. & 2nd. nn)
!         h01i = im.  part of q-parallel inter-plane hoppings (1st. & 2nd. nn)
!      
! set values of h00r, and h00i, h01r, and h01i equal to zero 

  h00r = 0.0d0
  h00i = 0.0d0         
  h01r = 0.0d0
  h01i = 0.0d0
  h02r = 0.d0
  h02i = 0.d0

  h00r(1,1) = s0(mm)
  forall (i=2:4) h00r(i,i) = p0(mm)
  forall (i=5:9) h00r(i,i) = d0(mm)

! scalar product k.r0(i); i=1,(n0/2) 


  n0 = n01 + n02

  do i=1,n0

! cosine and sine of k.r0(i); i=1,n0s2

   spkr0 = pi*( ax*r0(i,1) + ay*r0(i,2) )
   ckr0  = cos(spkr0)
   skr0  = sin(spkr0)

! calculation of h00r and h00i

   if (i>n01) then
    inn=2
   else
    inn=1
   end if
        
   g1=sss(mm,inn)
   g2=sps(mm,inn)
   g3=pps(mm,inn)
   g4=ppp(mm,inn)
   g5=sds(mm,inn)
   g6=pds(mm,inn)
   g7=pdp(mm,inn)
   g8=dds(mm,inn)
   g9=ddp(mm,inn)
   g10=ddd(mm,inn)
   w = c0(i,:)
   call eint(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,w,bp)
   h00r = h00r + ckr0*bp
   h00i = h00i + skr0*bp

  end do
 
! scalar product (k.r1(i)); i=1,n1

  do i=1,n1
 
! cosine(k.r1(i)) and sine(k.r1(i)); i=1,n1 

   spkr1 = pi*( ax*r1(i,1) + ay*r1(i,2) )
   ckr1  = cos(spkr1) 
   skr1  = dsin(spkr1)

! calculation of h01r and h01i
        
   g1=sss(mm,1)
   g2=sps(mm,1)
   g3=pps(mm,1)
   g4=ppp(mm,1)
   g5=sds(mm,1)
   g6=pds(mm,1)
   g7=pdp(mm,1)
   g8=dds(mm,1)
   g9=ddp(mm,1)
   g10=ddd(mm,1)
   w = c1(i,:)
   call eint(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,w,bp)
   h01r = h01r + ckr1*bp
   h01i = h01i + skr1*bp

  end do 

! h02: 
! scalar product (k.r2(i)); i=1,n2

  do i=1,n2
 
! cosine(k.r2(i)) and sine(k.r2(i)); i=1,n2 

   spkr2 = pi*( ax*r2(i,1) + ay*r2(i,2) )
   ckr2  = cos(spkr2) 
   skr2  = sin(spkr2)

! calculation of h02r and h02i
        
   g1=sss(mm,2)
   g2=sps(mm,2)
   g3=pps(mm,2)
   g4=ppp(mm,2)
   g5=sds(mm,2)
   g6=pds(mm,2)
   g7=pdp(mm,2)
   g8=dds(mm,2)
   g9=ddp(mm,2)
   g10=ddd(mm,2)
   w = c2(i,:)
   call eint(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,w,bp)
   h02r = h02r + ckr2*bp
   h02i = h02i + skr2*bp

  end do 

  return 

end subroutine helphfcc 
