SUBROUTINE green(wc,gout_Bi,gout_ad)

 USE f90_kind
 USE constants_u
 USE lattice
 USE the_hamiltonian
 USE the_LS_p_matrix
 USE the_LS_matrix
 USE splitting
 USE lambda_SOC

 COMPLEX(double) :: wc
 COMPLEX(double), DIMENSION(N_Bi,N_Bi,18,18), INTENT(OUT) :: gout_Bi
 COMPLEX(double), DIMENSION(Nadatoms,Nadatoms,18,18), INTENT(OUT) :: gout_ad

 COMPLEX(double), DIMENSION(dimG,dimG) :: g0i, gaux

 COMPLEX(double), DIMENSION(2*dimH_p3,2*dimH_p3) :: g0_l,gl_l,gl_r,gl_b

 REAL(double) :: kdotb
 COMPLEX(double) :: eikb,emikb

 INTEGER :: is,i,j,mu,mu0,muf,nu0,nuf,ineighb

 INTEGER :: dimdecim,itmax
 REAL(double) :: cri

 ! Decim
 cri = 1d-5
 itmax = 20

 g0_l = CMPLX(0d0,0d0,double)

 g0_l(1:dimH_p3,1:dimH_p3) = H00_pristine
 g0_l(dimH_p3+1:2*dimH_p3,dimH_p3+1:2*dimH_p3) = g0_l(1:dimH_p3,1:dimH_p3)

 ! setting up SOC in Bi leads
 DO is=1,N_Bi
 mu0 = (is-1)*9+2; muf = mu0 + 2
 nu0 = (is-1)*9+2; nuf = nu0 + 2
 g0_l(mu0:muf,nu0:nuf) = g0_l(mu0:muf,nu0:nuf) + lambda_Bi_p*LS_p(1:3,1:3)
 mu0 = (is-1)*9+5; muf = mu0 + 4
 nu0 = (is-1)*9+5; nuf = nu0 + 4
 g0_l(mu0:muf,nu0:nuf) = g0_l(mu0:muf,nu0:nuf) + lambda_Bi_d*LS(1,5:9,5:9)

 mu0 = (is-1)*9+2 + dimH_p3; muf = mu0 + 2
 nu0 = (is-1)*9+2 + dimH_p3; nuf = nu0 + 2
 g0_l(mu0:muf,nu0:nuf) = g0_l(mu0:muf,nu0:nuf) + lambda_Bi_p*LS_p(4:6,4:6)
 mu0 = (is-1)*9+5 + dimH_p3; muf = mu0 + 4
 nu0 = (is-1)*9+5 + dimH_p3; nuf = nu0 + 4
 g0_l(mu0:muf,nu0:nuf) = g0_l(mu0:muf,nu0:nuf) + lambda_Bi_d*LS(1,14:18,14:18)

 mu0 = (is-1)*9+2 + dimH_p3; muf = mu0 + 2
 nu0 = (is-1)*9+2 ; nuf = nu0 + 2
 g0_l(mu0:muf,nu0:nuf) = g0_l(mu0:muf,nu0:nuf) + lambda_Bi_p*LS_p(4:6,1:3)
 mu0 = (is-1)*9+5 + dimH_p3; muf = mu0 + 4
 nu0 = (is-1)*9+5 ; nuf = nu0 + 4
 g0_l(mu0:muf,nu0:nuf) = g0_l(mu0:muf,nu0:nuf) + lambda_Bi_d*LS(1,14:18,5:9)

 mu0 = (is-1)*9+2 ; muf = mu0 + 2
 nu0 = (is-1)*9+2 + dimH_p3; nuf = nu0 + 2
 g0_l(mu0:muf,nu0:nuf) = g0_l(mu0:muf,nu0:nuf) + lambda_Bi_p*LS_p(1:3,4:6)
 mu0 = (is-1)*9+5 ; muf = mu0 + 4
 nu0 = (is-1)*9+5 + dimH_p3; nuf = nu0 + 4
 g0_l(mu0:muf,nu0:nuf) = g0_l(mu0:muf,nu0:nuf) + lambda_Bi_d*LS(1,5:9,14:18)
 END DO

 g0_l = wc*mum_l - g0_l

 dimdecim = 2*dimH_p3
 CALL decim(g0_l,tl,tld,gl_l,gl_r,gl_b,cri,dimdecim,itmax)

 CALL invers(gl_l,dimdecim)
 CALL invers(gl_r,dimdecim)

 ! Central region (containing the adatom)

 g0i = CMPLX(0.d0,0.d0,double)

 g0i(1:dimH,1:dimH) = H00_central
 g0i(dimH+1:dimG,dimH+1:dimG) = g0i(1:dimH,1:dimH)

 ! Setting up SOC in the central region and
 ! spin splitting in the adatom
 ! Atom 1: adatom
 ! Atoms 2 to N_Bi+1: Bi
 ! last 2 atoms: H
 mu0 = Norb_ad - 4; muf = Norb_ad
 nu0 = Norb_ad - 4; nuf = Norb_ad
 g0i(mu0:muf,nu0:nuf) = g0i(mu0:muf,nu0:nuf) + lambda_ad*LS(1,5:9,5:9)
 mu0 = dimH + Norb_ad - 4; muf = dimH + Norb_ad
 nu0 = dimH + Norb_ad - 4; nuf = dimH + Norb_ad
 g0i(mu0:muf,nu0:nuf) = g0i(mu0:muf,nu0:nuf) + lambda_ad*LS(1,14:18,14:18)
 mu0 = dimH + Norb_ad - 4; muf = dimH + Norb_ad
 nu0 = Norb_ad - 4; nuf = Norb_ad
 g0i(mu0:muf,nu0:nuf) = g0i(mu0:muf,nu0:nuf) + lambda_ad*LS(1,14:18,5:9)
 mu0 = Norb_ad - 4; muf = Norb_ad
 nu0 = dimH + Norb_ad - 4; nuf = dimH + Norb_ad
 g0i(mu0:muf,nu0:nuf) = g0i(mu0:muf,nu0:nuf) + lambda_ad*LS(1,5:9,14:18)

 DO is=1,N_Bi
 mu0 = Nadatoms*Norb_ad + (is-1)*9+2; muf = mu0 + 2
 nu0 = Nadatoms*Norb_ad + (is-1)*9+2; nuf = nu0 + 2
 g0i(mu0:muf,nu0:nuf) = g0i(mu0:muf,nu0:nuf) + lambda_Bi_p*LS_p(1:3,1:3)
 mu0 = Nadatoms*Norb_ad + (is-1)*9+5; muf = mu0 + 4
 nu0 = Nadatoms*Norb_ad + (is-1)*9+5; nuf = nu0 + 4
 g0i(mu0:muf,nu0:nuf) = g0i(mu0:muf,nu0:nuf) + lambda_Bi_d*LS(1,5:9,5:9)

 mu0 = Nadatoms*Norb_ad + (is-1)*9+2 + dimH; muf = mu0 + 2
 nu0 = Nadatoms*Norb_ad + (is-1)*9+2 + dimH; nuf = nu0 + 2
 g0i(mu0:muf,nu0:nuf) = g0i(mu0:muf,nu0:nuf) + lambda_Bi_p*LS_p(4:6,4:6)
 mu0 = Nadatoms*Norb_ad + (is-1)*9+5 + dimH; muf = mu0 + 4
 nu0 = Nadatoms*Norb_ad + (is-1)*9+5 + dimH; nuf = nu0 + 4
 g0i(mu0:muf,nu0:nuf) = g0i(mu0:muf,nu0:nuf) + lambda_Bi_d*LS(1,14:18,14:18)

 mu0 = Nadatoms*Norb_ad + (is-1)*9+2 + dimH; muf = mu0 + 2
 nu0 = Nadatoms*Norb_ad + (is-1)*9+2 ; nuf = nu0 + 2
 g0i(mu0:muf,nu0:nuf) = g0i(mu0:muf,nu0:nuf) + lambda_Bi_p*LS_p(4:6,1:3)
 mu0 = Nadatoms*Norb_ad + (is-1)*9+5 + dimH; muf = mu0 + 4
 nu0 = Nadatoms*Norb_ad + (is-1)*9+5 ; nuf = nu0 + 4
 g0i(mu0:muf,nu0:nuf) = g0i(mu0:muf,nu0:nuf) + lambda_Bi_d*LS(1,14:18,5:9)

 mu0 = Nadatoms*Norb_ad + (is-1)*9+2 ; muf = mu0 + 2
 nu0 = Nadatoms*Norb_ad + (is-1)*9+2 + dimH; nuf = nu0 + 2
 g0i(mu0:muf,nu0:nuf) = g0i(mu0:muf,nu0:nuf) + lambda_Bi_p*LS_p(1:3,4:6)
 mu0 = Nadatoms*Norb_ad + (is-1)*9+5 ; muf = mu0 + 4
 nu0 = Nadatoms*Norb_ad + (is-1)*9+5 + dimH; nuf = nu0 + 4
 g0i(mu0:muf,nu0:nuf) = g0i(mu0:muf,nu0:nuf) + lambda_Bi_d*LS(1,5:9,14:18)
 END DO

 DO mu=Norb_ad-4,Norb_ad
 g0i(mu,mu) = g0i(mu,mu) + CMPLX(e0d(1) - hdel_A(1) - hw0,0d0,double)
 g0i(mu+dimH,mu+dimH) = g0i(mu+dimH,mu+dimH) + &
                             CMPLX(e0d(1) + hdel_A(1) + hw0,0d0,double)
 END DO

 g0i = wc*mum - g0i


 gaux = g0i - MATMUL(Hcll,MATMUL(gl_r,Hllc)) - MATMUL(Hcrl,MATMUL(gl_l,Hrlc))
 
 CALL invers(gaux,dimG)


 DO i=1,N_Bi
 DO j=1,N_Bi
 mu0 = Norb_ad*Nadatoms + (i-1)*Norb_Bi + 1; muf = mu0 + Norb_Bi - 1
 nu0 = Norb_ad*Nadatoms + (j-1)*Norb_Bi + 1; nuf = nu0 + Norb_Bi - 1
 gout_Bi(i,j,1:Norb_bi,1:Norb_Bi) = gaux(mu0:muf,nu0:nuf)
 gout_Bi(i,j,1:Norb_Bi,Norb_Bi+1:2*Norb_Bi) = gaux(mu0:muf,dimH+nu0:dimH+nuf)
 gout_Bi(i,j,Norb_Bi+1:2*Norb_Bi,1:Norb_Bi) = gaux(dimH+mu0:dimH+muf,nu0:nuf)
 gout_Bi(i,j,Norb_Bi+1:2*Norb_Bi,Norb_Bi+1:2*Norb_Bi) = &
                                   gaux(dimH+mu0:dimH+muf,dimH+nu0:dimH+nuf)
 END DO
 END DO


 DO i=1,Nadatoms
 DO j=1,Nadatoms
 mu0 = (i-1)*Norb_ad + 1; muf = mu0 + Norb_ad - 1
 nu0 = (j-1)*Norb_ad + 1; nuf = nu0 + Norb_ad - 1
 gout_ad(i,j,1:Norb_ad,1:Norb_ad) = gaux(mu0:muf,nu0:nuf)
 gout_ad(i,j,1:Norb_ad,Norb_ad+1:2*Norb_ad) = gaux(mu0:muf,dimH+nu0:dimH+nuf)
 gout_ad(i,j,Norb_ad+1:2*Norb_ad,1:Norb_ad) = gaux(dimH+mu0:dimH+muf,nu0:nuf)
 gout_ad(i,j,Norb_ad+1:2*Norb_ad,Norb_ad+1:2*Norb_ad) = &
                                    gaux(dimH+mu0:dimH+muf,dimH+nu0:dimH+nuf)
 END DO
 END DO
 RETURN

END SUBROUTINE green

