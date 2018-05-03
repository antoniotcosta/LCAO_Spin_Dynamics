SUBROUTINE the_L_matrix(L_x,L_y,L_z)

 USE f90_kind

 COMPLEX(double), DIMENSION(9,9), INTENT(OUT) :: L_x,L_y,L_z

 COMPLEX(double), DIMENSION(18,18) :: LS 
 COMPLEX(double), DIMENSION(9,9) :: hL_p,hL_m

 REAL(double) :: theta,phi
 COMPLEX(double) :: zi 

 zi = CMPLX(0.d0,1.d0,double)

 theta = 0.d0
 phi = 0.d0

 CALL LS_matrix(theta,phi,LS)

 L_z = 2.d0*LS(1:9,1:9)

 hL_m = LS(1:9,10:18)
 hL_p = LS(10:18,1:9)

 L_x = hL_p + hL_m
 L_y = -zi*( hL_p - hL_m )

 RETURN

END SUBROUTINE the_L_matrix

