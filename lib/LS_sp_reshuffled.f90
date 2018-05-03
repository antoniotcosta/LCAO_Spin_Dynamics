SUBROUTINE LS_matrix_p(theta,phi)

 USE f90_kind
 USE the_LS_p_matrix

 REAL(double), INTENT(IN) :: theta,phi
 COMPLEX(double) :: z0,zi,z1

 COMPLEX(double), DIMENSION(2,2) :: Sx,Sy,Sz
 COMPLEX(double), DIMENSION(3,3) :: Lx,Ly,Lz

 COMPLEX(double), DIMENSION(6,6) :: LS_tmp

 INTEGER, DIMENSION(3) :: inew
 INTEGER :: i,j

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

  LS_p = z0 

  LS_p(1:3,1:3) = Sx(1,1)*Lx + Sy(1,1)*Ly + Sz(1,1)*Lz
  LS_p(4:6,4:6) = Sx(2,2)*Lx + Sy(2,2)*Ly + Sz(2,2)*Lz
  LS_p(1:3,4:6) = Sx(1,2)*Lx + Sy(1,2)*Ly + Sz(1,2)*Lz
  LS_p(4:6,1:3) = Sx(2,1)*Lx + Sy(2,1)*Ly + Sz(2,1)*Lz

  ! Reshuffling

  inew = [ 3, 1, 2 ]

  DO i=1,3
  DO j=1,3
  LS_tmp(  i,  j) = LS_p(inew(i),inew(j)) 
  LS_tmp(i+3,  j) = LS_p(inew(i)+3,inew(j)) 
  LS_tmp(i  ,j+3) = LS_p(inew(i),inew(j)+3) 
  LS_tmp(i+3,j+3) = LS_p(inew(i)+3,inew(j)+3) 
  END DO
  END DO

  LS_p = LS_tmp

  
 RETURN

END SUBROUTINE LS_matrix_p
