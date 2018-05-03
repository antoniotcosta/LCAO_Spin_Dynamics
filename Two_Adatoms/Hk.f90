MODULE H_of_k

CONTAINS
SUBROUTINE Hk(qv,gout)

 USE f90_kind
 USE constants_u
 USE lattice
 USE the_hamiltonian
 USE the_LS_p_matrix
 USE the_LS_matrix
 USE splitting
 USE lambda_SOC

 REAL(double), DIMENSION(2), INTENT(IN) :: qv
 COMPLEX(double), DIMENSION(:,:), INTENT(OUT) :: gout

 REAL(double) :: kdotb
 COMPLEX(double) :: eikb,emikb

 INTEGER :: is,i,j,mu,mu0,muf,nu0,nuf,ineighb

 gout = CMPLX(0.d0,0.d0,double)

 gout(1:dimH,1:dimH) = CMPLX(H00,0.d0,double)

 DO ineighb=1,numbneighb

 kdotb = 2.d0*pi*DOT_PRODUCT(qv,dij(ineighb,1:2))/dij(1,1)

 eikb = CMPLX( COS(kdotb), SIN(kdotb), double)
 emikb = CMPLX( COS(kdotb), -SIN(kdotb), double)


 gout(1:dimH,1:dimH) = gout(1:dimH,1:dimH) + & 
                      CMPLX(Hij(ineighb,:,:),0.d0,double)*eikb + &
                      CMPLX(TRANSPOSE(Hij(ineighb,:,:)),0.d0,double)*emikb
 END DO

 gout(dimH+1:dimG,dimH+1:dimG) = gout(1:dimH,1:dimH)

 ! Atom 1: adatom
 ! Atoms 2 to N_Bi+1: Bi
 ! last 2 atoms: H
 mu0 = Norb_ad - 4; muf = Norb_ad
 nu0 = Norb_ad - 4; nuf = Norb_ad
 gout(mu0:muf,nu0:nuf) = gout(mu0:muf,nu0:nuf) + lambda_ad*LS(5:9,5:9)
 mu0 = dimH + Norb_ad - 4; muf = dimH + Norb_ad
 nu0 = dimH + Norb_ad - 4; nuf = dimH + Norb_ad
 gout(mu0:muf,nu0:nuf) = gout(mu0:muf,nu0:nuf) + lambda_ad*LS(14:18,14:18)
 mu0 = dimH + Norb_ad - 4; muf = dimH + Norb_ad
 nu0 = Norb_ad - 4; nuf = Norb_ad
 gout(mu0:muf,nu0:nuf) = gout(mu0:muf,nu0:nuf) + lambda_ad*LS(14:18,5:9)
 mu0 = Norb_ad - 4; muf = Norb_ad
 nu0 = dimH + Norb_ad - 4; nuf = dimH + Norb_ad
 gout(mu0:muf,nu0:nuf) = gout(mu0:muf,nu0:nuf) + lambda_ad*LS(5:9,14:18)

 DO is=1,N_Bi
 mu0 = Nadatoms*Norb_ad + (is-1)*9+2; muf = mu0 + 2
 nu0 = Nadatoms*Norb_ad + (is-1)*9+2; nuf = nu0 + 2
 gout(mu0:muf,nu0:nuf) = gout(mu0:muf,nu0:nuf) + lambda_Bi_p*LS_p(2:4,2:4)
 mu0 = Nadatoms*Norb_ad + (is-1)*9+5; muf = mu0 + 4
 nu0 = Nadatoms*Norb_ad + (is-1)*9+5; nuf = nu0 + 4
 gout(mu0:muf,nu0:nuf) = gout(mu0:muf,nu0:nuf) + lambda_Bi_d*LS(5:9,5:9)

 mu0 = Nadatoms*Norb_ad + (is-1)*9+2 + dimH; muf = mu0 + 2
 nu0 = Nadatoms*Norb_ad + (is-1)*9+2 + dimH; nuf = nu0 + 2
 gout(mu0:muf,nu0:nuf) = gout(mu0:muf,nu0:nuf) + lambda_Bi_p*LS_p(6:8,6:8)
 mu0 = Nadatoms*Norb_ad + (is-1)*9+5 + dimH; muf = mu0 + 4
 nu0 = Nadatoms*Norb_ad + (is-1)*9+5 + dimH; nuf = nu0 + 4
 gout(mu0:muf,nu0:nuf) = gout(mu0:muf,nu0:nuf) + lambda_Bi_d*LS(14:18,14:18)

 mu0 = Nadatoms*Norb_ad + (is-1)*9+2 + dimH; muf = mu0 + 2
 nu0 = Nadatoms*Norb_ad + (is-1)*9+2 ; nuf = nu0 + 2
 gout(mu0:muf,nu0:nuf) = gout(mu0:muf,nu0:nuf) + lambda_Bi_p*LS_p(6:8,2:4)
 mu0 = Nadatoms*Norb_ad + (is-1)*9+5 + dimH; muf = mu0 + 4
 nu0 = Nadatoms*Norb_ad + (is-1)*9+5 ; nuf = nu0 + 4
 gout(mu0:muf,nu0:nuf) = gout(mu0:muf,nu0:nuf) + lambda_Bi_d*LS(14:18,5:9)

 mu0 = Nadatoms*Norb_ad + (is-1)*9+2 ; muf = mu0 + 2
 nu0 = Nadatoms*Norb_ad + (is-1)*9+2 + dimH; nuf = nu0 + 2
 gout(mu0:muf,nu0:nuf) = gout(mu0:muf,nu0:nuf) + lambda_Bi_p*LS_p(2:4,6:8)
 mu0 = Nadatoms*Norb_ad + (is-1)*9+5 ; muf = mu0 + 4
 nu0 = Nadatoms*Norb_ad + (is-1)*9+5 + dimH; nuf = nu0 + 4
 gout(mu0:muf,nu0:nuf) = gout(mu0:muf,nu0:nuf) + lambda_Bi_d*LS(5:9,14:18)
 END DO

 DO mu=Norb_ad-4,Norb_ad
 gout(mu,mu) = gout(mu,mu) + CMPLX(e0d(1) - hdel_A(1) - hw0,0d0,double)
 gout(mu+dimH,mu+dimH) = gout(mu+dimH,mu+dimH) + &
                             CMPLX(e0d(1) + hdel_A(1) + hw0,0d0,double)
 END DO



 RETURN

END SUBROUTINE Hk

END MODULE H_of_k
