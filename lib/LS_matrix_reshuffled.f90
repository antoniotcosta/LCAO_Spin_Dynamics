SUBROUTINE LS_matrix(theta,phi,LS)

 USE f90_kind
 USE constants_u

 COMPLEX(double), DIMENSION(18,18), INTENT(OUT) :: LS
 COMPLEX(double), DIMENSION(18,18) :: LS_tmp
 COMPLEX(double), DIMENSION(5,5) :: L_x,L_y,L_z

 REAL(double), INTENT(IN) :: theta,phi

 INTEGER :: i,j
 INTEGER, DIMENSION(9) :: inew



  ! The Spin-Orbit matrix

  LS = zero

  ! The d-block
  LS( 5, 6) =  0.5d0*zi*SIN(theta)*SIN(phi) 
  LS( 5, 7) = -0.5d0*zi*SIN(theta)*COS(phi)
  LS( 5, 8) =  zi*COS(theta)

  LS( 6, 5) =  CONJG(LS(5,6))
  LS( 6, 7) =  0.5d0*zi*COS(theta)
  LS( 6, 8) = -0.5d0*zi*SIN(theta)*COS(phi)
  LS( 6, 9) = -0.5d0*zi*sq3*SIN(theta)*COS(phi)
  
  LS( 7, 5) =  CONJG(LS(5,7))
  LS( 7, 6) =  CONJG(LS(6,7))
  LS( 7, 8) = -0.5d0*zi*SIN(theta)*SIN(phi)
  LS( 7, 9) =  0.5d0*zi*sq3*SIN(theta)*SIN(phi)
  
  LS( 8, 5) =  CONJG(LS(5,8))
  LS( 8, 6) =  CONJG(LS(6,8))
  LS( 8, 7) =  CONJG(LS(7,8))
  
  LS( 9, 6) =  CONJG(LS(6,9))
  LS( 9, 7) =  CONJG(LS(7,9))


  LS( 5,15) =  0.5d0*(COS(phi)+zi*COS(theta)*SIN(phi))
  LS( 5,16) =  0.5d0*(SIN(phi)-zi*COS(theta)*COS(phi))
  LS( 5,17) = -zi*SIN(theta)

  LS( 6,14) = -0.5d0*(COS(phi)+zi*COS(theta)*SIN(phi))
  LS( 6,16) = -0.5d0*zi*SIN(theta)
  LS( 6,17) =  0.5d0*(SIN(phi)-zi*COS(theta)*COS(phi))
  LS( 6,18) =  0.5d0*sq3*(SIN(phi)-zi*COS(theta)*COS(phi))

  LS( 7,14) =  0.5D0*(-SIN(phi)+zi*COS(theta)*COS(phi))
  LS( 7,15) =  0.5d0*zi*SIN(theta)
  LS( 7,17) = -0.5d0*(COS(phi)+zi*COS(theta)*SIN(phi))
  LS( 7,18) =  0.5d0*sq3*(COS(phi)+zi*COS(theta)*SIN(phi))

  LS( 8,14) =  zi*SIN(theta)
  LS( 8,15) =  0.5d0*(-SIN(phi)+zi*COS(theta)*COS(phi))
  LS( 8,16) =  0.5d0*(COS(phi)+zi*COS(theta)*SIN(phi))

  LS( 9,15) =  0.5d0*sq3*(-SIN(phi)+zi*COS(theta)*COS(phi))
  LS( 9,16) = -0.5d0*sq3*(COS(phi)+zi*COS(theta)*SIN(phi))

  LS(14:18,14:18) = -LS(5:9,5:9)

  LS(14:18,5:9) = TRANSPOSE(CONJG(LS(5:9,14:18)))

  ! Reshuffling

  inew = [ 1, 4, 2, 3, 9, 7, 6, 8, 5 ]

  DO i=1,9
  DO j=1,9
  LS_tmp(i,j) = LS(inew(i),inew(j)) 
  LS_tmp(i+9,j) = LS(inew(i)+9,inew(j)) 
  LS_tmp(i,j+9) = LS(inew(i),inew(j)+9) 
  LS_tmp(i+9,j+9) = LS(inew(i)+9,inew(j)+9) 
  END DO
  END DO

  LS = LS_tmp

  L_x = LS(14:18, 5: 9) + LS( 5: 9,14:18) 
  L_y = -zi*( LS(14:18, 5: 9) -  LS( 5: 9,14:18) )
  L_z = 2.d0*LS( 5: 9, 5: 9) 

  RETURN

END SUBROUTINE LS_matrix
