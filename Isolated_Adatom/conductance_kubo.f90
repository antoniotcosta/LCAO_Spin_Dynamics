SUBROUTINE conductance_kubo(G00,G11,G01,G10,H01,H10,Gamma_Tr,Gamma_per_atom)
! Implements expression (7) from PRB 55, 14378 (1997).

 USE f90_kind
 USE lattice, ONLY : dimG,dimH_p3,N_uc

 COMPLEX(double), DIMENSION(dimG,dimG), INTENT(IN) :: G00
 COMPLEX(double), DIMENSION(2*dimH_p3,2*dimH_p3), INTENT(IN) :: G11
 COMPLEX(double), DIMENSION(dimG,2*dimH_p3), INTENT(IN) :: G01,H01
 COMPLEX(double), DIMENSION(2*dimH_p3,dimG), INTENT(IN) :: G10,H10
 REAL(double), INTENT(OUT) :: Gamma_Tr
 REAL(double), DIMENSION(N_uc-6), INTENT(OUT) :: Gamma_per_atom

 COMPLEX(double), DIMENSION(dimG,dimG) :: Gtilde00
 COMPLEX(double), DIMENSION(2*dimH_p3,2*dimH_p3) :: Gtilde11
 COMPLEX(double), DIMENSION(2*dimH_p3,dimG) :: Gtilde10
 COMPLEX(double), DIMENSION(dimG,dimG) :: Gamma_M

 COMPLEX(double) :: zi

 INTEGER :: mu

 zi = CMPLX(0d0,1d0,double)

 Gtilde00 = .5d0*zi*( G00 - TRANSPOSE(CONJG(G00)) )
 Gtilde11 = .5d0*zi*( G11 - TRANSPOSE(CONJG(G11)) )
 Gtilde10 = .5d0*zi*( G10 - TRANSPOSE(CONJG(G01)) )

 Gamma_M = MATMUL(MATMUL(Gtilde00,H01),MATMUL(Gtilde11,H10)) - &
           MATMUL(MATMUL(H01,Gtilde10),MATMUL(H01,Gtilde10)) 

 Gamma_Tr = 0d0
 DO mu=1,dimG
 Gamma_Tr = Gamma_Tr + 4d0*REAL(Gamma_M(mu,mu))
 END DO

 Gamma_per_atom = 0d0
 DO i=1,N_uc-6
 DO mu=1,9
 mup = (i-1)*9+mu
 Gamma_per_atom(i) = Gamma_per_atom(i) + 4d0*REAL(Gamma_M(mup,mup))
 END DO
 DO mu=1,9
 mup = (i-1)*9+mu+dimG
 Gamma_per_atom(i) = Gamma_per_atom(i) + 4d0*REAL(Gamma_M(mup,mup))
 END DO
 END DO

 RETURN

END SUBROUTINE conductance_kubo
