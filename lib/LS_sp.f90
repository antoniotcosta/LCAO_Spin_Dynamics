SUBROUTINE LS_matrix_sp(theta,phi)

 USE f90_kind
 USE the_LS_sp_matrix

 REAL(double), INTENT(IN) :: theta,phi
 COMPLEX(double) :: z0,zi,z1

 COMPLEX(double), DIMENSION(2,2) :: Sx,Sy,Sz
 COMPLEX(double), DIMENSION(3,3) :: Lx,Ly,Lz

 z0 = CMPLX(0d0,0d0,double)
 zi = CMPLX(0d0,1d0,double)
 z1 = CMPLX(1d0,0d0,double)

 ! S matrices in arbitrary coordinate system

  Sx(1,1) = .5d0*CMPLX(SIN(theta)*COS(phi),0.d0,double) 
  Sx(2,2) = -Sx(1,1)
  Sx(1,2) = .5d0*CMPLX(COS(theta)*COS(phi),SIN(phi),double)
  Sx(2,1) = CONJG(Sx(1,2))

  Sy(1,1) = .5d0*CMPLX(SIN(theta)*SIN(phi),0.d0,double)
  Sy(2,2) = -Sy(1,1)
  Sy(1,2) = .5d0*CMPLX(COS(theta)*SIN(phi),-COS(phi),double)
  Sy(2,1) = CONJG(Sy(1,2))

  Sz(1,1) = .5d0*CMPLX(COS(theta),0.d0,double)
  Sz(2,2) = -Sz(1,1)
  Sz(1,2) = .5d0*CMPLX(-SIN(theta),0.d0,double)
  Sz(2,1) = Sz(1,2)

 ! L matrices (l=1) in standard coordinate system

  Lx(1,:) = z0
  Lx(2,:) = [z0, z0, -zi]
  Lx(3,:) = [z0, zi, z0]

  Ly(1,:) = [z0, z0, zi]
  Ly(2,:) = z0
  Ly(3,:) = [-zi, z0, z0]

  Lz(1,:) = [z0, -zi, z0]
  Lz(2,:) = [zi, z0, z0]
  Lz(3,:) = z0

  LS_sp = z0 

  LS_sp(2:4,2:4) = Sx(1,1)*Lx + Sy(1,1)*Ly + Sz(1,1)*Lz
  LS_sp(6:8,6:8) = Sx(2,2)*Lx + Sy(2,2)*Ly + Sz(2,2)*Lz
  LS_sp(2:4,6:8) = Sx(1,2)*Lx + Sy(1,2)*Ly + Sz(1,2)*Lz
  LS_sp(6:8,2:4) = Sx(2,1)*Lx + Sy(2,1)*Ly + Sz(2,1)*Lz

 RETURN

END SUBROUTINE LS_matrix_sp
