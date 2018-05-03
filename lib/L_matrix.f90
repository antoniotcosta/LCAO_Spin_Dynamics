SUBROUTINE L_matrix(Lx,Ly,Lz)

 USE f90_kind
 USE constants_u

 COMPLEX(double), DIMENSION(9,9), INTENT(OUT) :: Lx,Ly,Lz

 COMPLEX(double), DIMENSION(9,9) :: Lp,Lm

 Lz = zero

 Lz(2,3) = zi
 Lz(3,2) = -zi

 Lz(5,8) = -2.d0*zi
 Lz(8,5) = 2.d0*zi
 Lz(6,7) = -zi
 Lz(7,6) = zi

 Lp = zero
 Lm = zero

 Lp(4,2) = zum
 Lp(4,3) = -zi
 Lp(2,4) = -zum
 Lp(4,3) = zi

 Lp(5,6) = zum
 Lp(5,7) = -zi

 Lp(6,5) = -zum
 Lp(6,8) = -zi
 Lp(6,9) = -sq3*zi

 Lp(7,5) = zi
 Lp(7,8) = -zum
 Lp(7,9) = sq3

 Lp(8,6) = zi
 Lp(8,7) = zum

 Lp(9,6) = -sq3*zi
 Lp(9,7) = sq3 

 Lm = TRANSPOSE(CONJG(Lp))

 Lx = 0.5d0*(Lp+Lm)

 Ly = -0.5d0*zi*(Lp-Lm)

 RETURN

END SUBROUTINE L_matrix

