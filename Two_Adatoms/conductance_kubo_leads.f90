SUBROUTINE conductance_kubo(G00,G11,G01,G10,H01,H10,Gamma_Tr)
! Implements expression (7) from PRB 55, 14378 (1997).

 USE f90_kind
 USE lattice, ONLY : dimG,dimH_p3

 REAL(double), INTENT(OUT) :: Gamma_Tr
 COMPLEX(double), DIMENSION(2*dimH_p3,2*dimH_p3), INTENT(IN) :: G00,G11,G01,G10,H01,H10

 COMPLEX(double), DIMENSION(2*dimH_p3,2*dimH_p3) :: Gtilde11,Gtilde00,Gtilde10,Gamma_M

 COMPLEX(double) :: zi

 INTEGER :: mu

 zi = CMPLX(0d0,1d0,double)

 Gtilde00 = .5d0*zi*( G00 - TRANSPOSE(CONJG(G00)) )
 Gtilde11 = .5d0*zi*( G11 - TRANSPOSE(CONJG(G11)) )
 Gtilde10 = .5d0*zi*( G10 - TRANSPOSE(CONJG(G01)) )

 Gamma_M = MATMUL(MATMUL(Gtilde00,H01),MATMUL(Gtilde11,H10)) - &
           MATMUL(MATMUL(H01,Gtilde10),MATMUL(H01,Gtilde10)) 

 Gamma_Tr = 0d0
 DO mu=1,2*dimH_p3
 Gamma_Tr = Gamma_Tr + 4d0*REAL(Gamma_M(mu,mu))
 END DO

 RETURN

END SUBROUTINE conductance_kubo
