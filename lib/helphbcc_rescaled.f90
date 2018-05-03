SUBROUTINE helphbcc(ax,ay,az,h00r,h00i,h01r,h01i,mm)

 USE f90_kind
 USE constants_u
 USE lattice
 USE tight_binding

 INTEGER :: i,inn,mm,mml
 REAL(double), DIMENSION(9,9) :: h00r,h00i,h01r,h01i,bp
 REAL(double) :: w(3),ax,ay,az,spkr0,ckr0,skr0,spkr1,ckr1,skr1,spkr2,ckr2,skr2, &
                 g1,g2,g3,g4,g5,g6,g7,g8,g9,g10


!     input:
!         ax, ay, az = components of the wave vector k in units of 
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

  h00r(1,1) = s0(mm)
  FORALL (i=2:4) h00r(i,i) = p0(mm)
  FORALL (i=5:7) h00r(i,i) = d0t(mm)
  FORALL (i=8:9) h00r(i,i) = d0e(mm)             



  n0 = n01 + n02

  DO i=1,n0

     ! scalar product k.r0(i); i=1,(n0/2) 
     ! cosine and sine of k.r0(i); i=1,n0s2

     spkr0 = pi*( ax*r0(i,1) + ay*r0(i,2) + az*r0(i,3) )
     ckr0  = cos(spkr0)
     skr0  = sin(spkr0)

     ! calculation of h00r and h00i

     IF (i>n01) THEN
        inn=2
     ELSE
        inn=1
     END IF
        
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
     CALL eint(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,w,bp)
     h00r = h00r + ckr0*bp
     h00i = h00i + skr0*bp

   END DO
 

  DO i=1,n1
 
     ! scalar product (k.r1(i)); i=1,n1
     ! cosine(k.r1(i)) and sine(k.r1(i)); i=1,n1 

     spkr1 = pi*( ax*r1(i,1) + ay*r1(i,2) + az*r1(i,3) )
     ckr1  = cos(spkr1) 
     skr1  = dsin(spkr1)

     ! calculation of h01r and h01i
     mml = mm    
     IF (mm==2) mml = 3 ! we're rescaling the in-plane hopping,
     g1=sss(mml,1)      ! inter-plane hopping should be the same
     g2=sps(mml,1)
     g3=pps(mml,1)
     g4=ppp(mml,1)
     g5=sds(mml,1)
     g6=pds(mml,1)
     g7=pdp(mml,1)
     g8=dds(mml,1)
     g9=ddp(mml,1)
     g10=ddd(mml,1)
     w = c1(i,:)
     CALL eint(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,w,bp)
     h01r = h01r + ckr1*bp
     h01i = h01i + skr1*bp

  END DO 
 

  DO i=1,n2
 
     ! scalar product (k.r2(i)); i=1,n2
     ! cosine(k.r2(i)) and sine(k.r2(i)); i=1,n2 

     spkr2 = pi*( ax*r2(i,1) + ay*r2(i,2) + az*r2(i,3) )
     ckr2  = COS(spkr2) 
     skr2  = SIN(spkr2)

     ! calculation of h01r and h01i
        
     g1=sss(mml,2)
     g2=sps(mml,2)
     g3=pps(mml,2)
     g4=ppp(mml,2)
     g5=sds(mml,2)
     g6=pds(mml,2)
     g7=pdp(mml,2)
     g8=dds(mml,2)
     g9=ddp(mml,2)
     g10=ddd(mml,2)
     w = c2(i,:)
     CALL eint(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,w,bp)
     h01r = h01r + ckr2*bp
     h01i = h01i + skr2*bp

  END DO 
 
  RETURN 

END SUBROUTINE helphbcc 
