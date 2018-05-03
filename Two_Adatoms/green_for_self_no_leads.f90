! This is a test version with two adatom cells
! separated by a number of pristine cells, 
! without the leads.

SUBROUTINE green(wc,gout_ad)

 USE f90_kind
 USE constants_u
 USE lattice
 USE the_hamiltonian
 USE the_LS_p_matrix
 USE the_LS_matrix
 USE splitting
 USE lambda_SOC
 USE MPI_pars

 COMPLEX(double) :: wc
 COMPLEX(double), DIMENSION(Nadtotal,Nadtotal,18,18), INTENT(OUT) :: gout_ad

 COMPLEX(double), DIMENSION(dimG,dimG) :: g0i_A,g0i_B,gaux,g0AA,g0BB,gAA,gBB,gAB,gBA
 COMPLEX(double), DIMENSION(2*dimH_p3,2*dimH_p3) :: g0i_Pristine,g0Pristine
 COMPLEX(double), DIMENSION(2*dimH_p3,2*dimH_p3) :: g0Nm1Nm1,gNN
 COMPLEX(double), DIMENSION(dimG,2*dimH_p3) :: g0ANm1,gAN
 COMPLEX(double), DIMENSION(2*dimH_p3,dimG) :: g0Nm1A,gNA


 INTEGER :: is,i,j,mu,mu0,muf,nu0,nuf,ineighb,ispacer

 REAL(double) :: deltaG

 mum = CMPLX(0d0,0d0,double)
 DO i=1,dimG
 mum(i,i) = CMPLX(1d0,0d0,double) 
 END DO

 ! Central region (containing the adatoms)
 ! The central region consists of two adjacent
 ! cells with adatoms. The hopping between these 
 ! two adjacent cells is extracted from the DFT
 ! calculation for the periodic chain of adatoms.

 g0i_A = CMPLX(0.d0,0.d0,double)
 g0i_B = CMPLX(0.d0,0.d0,double)

 g0i_A(1:dimH,1:dimH) = H00_central
 g0i_A(dimH+1:dimG,dimH+1:dimG) = H00_central

 ! Setting up SOC in the central region and
 ! spin splitting in the adatom
 ! Atom 1: adatom
 ! Atoms 2 to N_Bi+1: Bi
 ! last 2 atoms: H
 mu0 = Norb_ad - 4; muf = Norb_ad
 nu0 = Norb_ad - 4; nuf = Norb_ad
 g0i_A(mu0:muf,nu0:nuf) = g0i_A(mu0:muf,nu0:nuf) + lambda_ad*LS(5:9,5:9)
 mu0 = dimH + Norb_ad - 4; muf = dimH + Norb_ad
 nu0 = dimH + Norb_ad - 4; nuf = dimH + Norb_ad
 g0i_A(mu0:muf,nu0:nuf) = g0i_A(mu0:muf,nu0:nuf) + lambda_ad*LS(14:18,14:18)
 mu0 = dimH + Norb_ad - 4; muf = dimH + Norb_ad
 nu0 = Norb_ad - 4; nuf = Norb_ad
 g0i_A(mu0:muf,nu0:nuf) = g0i_A(mu0:muf,nu0:nuf) + lambda_ad*LS(14:18,5:9)
 mu0 = Norb_ad - 4; muf = Norb_ad
 nu0 = dimH + Norb_ad - 4; nuf = dimH + Norb_ad
 g0i_A(mu0:muf,nu0:nuf) = g0i_A(mu0:muf,nu0:nuf) + lambda_ad*LS(5:9,14:18)

 DO is=1,N_Bi
 mu0 = Nadatoms*Norb_ad + (is-1)*9+2; muf = mu0 + 2
 nu0 = Nadatoms*Norb_ad + (is-1)*9+2; nuf = nu0 + 2
 g0i_A(mu0:muf,nu0:nuf) = g0i_A(mu0:muf,nu0:nuf) + lambda_Bi_p*LS_p(1:3,1:3)
 mu0 = Nadatoms*Norb_ad + (is-1)*9+5; muf = mu0 + 4
 nu0 = Nadatoms*Norb_ad + (is-1)*9+5; nuf = nu0 + 4
 g0i_A(mu0:muf,nu0:nuf) = g0i_A(mu0:muf,nu0:nuf) + lambda_Bi_d*LS(5:9,5:9)

 mu0 = Nadatoms*Norb_ad + (is-1)*9+2 + dimH; muf = mu0 + 2
 nu0 = Nadatoms*Norb_ad + (is-1)*9+2 + dimH; nuf = nu0 + 2
 g0i_A(mu0:muf,nu0:nuf) = g0i_A(mu0:muf,nu0:nuf) + lambda_Bi_p*LS_p(4:6,4:6)
 mu0 = Nadatoms*Norb_ad + (is-1)*9+5 + dimH; muf = mu0 + 4
 nu0 = Nadatoms*Norb_ad + (is-1)*9+5 + dimH; nuf = nu0 + 4
 g0i_A(mu0:muf,nu0:nuf) = g0i_A(mu0:muf,nu0:nuf) + lambda_Bi_d*LS(14:18,14:18)

 mu0 = Nadatoms*Norb_ad + (is-1)*9+2 + dimH; muf = mu0 + 2
 nu0 = Nadatoms*Norb_ad + (is-1)*9+2 ; nuf = nu0 + 2
 g0i_A(mu0:muf,nu0:nuf) = g0i_A(mu0:muf,nu0:nuf) + lambda_Bi_p*LS_p(4:6,1:3)
 mu0 = Nadatoms*Norb_ad + (is-1)*9+5 + dimH; muf = mu0 + 4
 nu0 = Nadatoms*Norb_ad + (is-1)*9+5 ; nuf = nu0 + 4
 g0i_A(mu0:muf,nu0:nuf) = g0i_A(mu0:muf,nu0:nuf) + lambda_Bi_d*LS(14:18,5:9)

 mu0 = Nadatoms*Norb_ad + (is-1)*9+2 ; muf = mu0 + 2
 nu0 = Nadatoms*Norb_ad + (is-1)*9+2 + dimH; nuf = nu0 + 2
 g0i_A(mu0:muf,nu0:nuf) = g0i_A(mu0:muf,nu0:nuf) + lambda_Bi_p*LS_p(1:3,4:6)
 mu0 = Nadatoms*Norb_ad + (is-1)*9+5 ; muf = mu0 + 4
 nu0 = Nadatoms*Norb_ad + (is-1)*9+5 + dimH; nuf = nu0 + 4
 g0i_A(mu0:muf,nu0:nuf) = g0i_A(mu0:muf,nu0:nuf) + lambda_Bi_d*LS(5:9,14:18)
 END DO

 g0i_B = g0i_A

 DO mu=Norb_ad-4,Norb_ad
 g0i_A(mu,mu) = g0i_A(mu,mu) + CMPLX(e0d(1)-hw0,0d0,double) + hdel_A(1,1,1) 
 g0i_A(mu+dimH,mu+dimH) = g0i_A(mu+dimH,mu+dimH) + &
                               CMPLX(e0d(1)+hw0,0d0,double) + hdel_A(1,2,2) 
 g0i_A(mu,mu+dimH) = g0i_A(mu,mu+dimH) + hdel_A(1,1,2) 
 g0i_A(mu+dimH,mu) = g0i_A(mu+dimH,mu) + hdel_A(1,2,1) 

 g0i_B(mu,mu) = g0i_B(mu,mu) + CMPLX(e0d(2)-hw0,0d0,double) + hdel_A(2,1,1) 
 g0i_B(mu+dimH,mu+dimH) = g0i_B(mu+dimH,mu+dimH) + &
                               CMPLX(e0d(2)+hw0,0d0,double) + hdel_A(2,2,2) 
 g0i_B(mu,mu+dimH) = g0i_B(mu,mu+dimH) + hdel_A(2,1,2) 
 g0i_B(mu+dimH,mu) = g0i_B(mu+dimH,mu) + hdel_A(2,2,1) 
 END DO

 g0i_A = wc*mum - g0i_A
 g0i_B = wc*mum - g0i_B

 ! Up to here we have to 3x1 flakes with one Fe adatom.
 ! The adatom is located at the left edge of the cell,
 ! thus we need to attach a pristine flake to to left (A)
 ! flake to make the A adatom "more similar" to the B adatom.
 ! Just to be sure, we also attach a pristine flake to the
 ! flake at the right.

 g0Pristine = CMPLX(0.d0,0.d0,double)
 g0Pristine(1:dimH_p3,1:dimH_p3) = H00_pristine
 g0Pristine(dimH_p3+1:2*dimH_p3,dimH_p3+1:2*dimH_p3) = H00_pristine

 DO is=1,N_Bi
 mu0 = (is-1)*9+2; muf = mu0 + 2
 nu0 = (is-1)*9+2; nuf = nu0 + 2
 g0Pristine(mu0:muf,nu0:nuf) = g0Pristine(mu0:muf,nu0:nuf) + &
                                 lambda_Bi_p*LS_p(1:3,1:3)
 mu0 = (is-1)*9+5; muf = mu0 + 4
 nu0 = (is-1)*9+5; nuf = nu0 + 4
 g0Pristine(mu0:muf,nu0:nuf) = g0Pristine(mu0:muf,nu0:nuf) + &
                                 lambda_Bi_d*LS(5:9,5:9)

 mu0 = (is-1)*9+2 + dimH_p3; muf = mu0 + 2
 nu0 = (is-1)*9+2 + dimH_p3; nuf = nu0 + 2
 g0Pristine(mu0:muf,nu0:nuf) = g0Pristine(mu0:muf,nu0:nuf) + &
                                 lambda_Bi_p*LS_p(4:6,4:6)
 mu0 = (is-1)*9+5 + dimH_p3; muf = mu0 + 4
 nu0 = (is-1)*9+5 + dimH_p3; nuf = nu0 + 4
 g0Pristine(mu0:muf,nu0:nuf) = g0Pristine(mu0:muf,nu0:nuf) + &
                                 lambda_Bi_d*LS(14:18,14:18)

 mu0 = (is-1)*9+2 + dimH_p3; muf = mu0 + 2
 nu0 = (is-1)*9+2 ; nuf = nu0 + 2
 g0Pristine(mu0:muf,nu0:nuf) = g0Pristine(mu0:muf,nu0:nuf) + &
                                 lambda_Bi_p*LS_p(4:6,1:3)
 mu0 = (is-1)*9+5 + dimH_p3; muf = mu0 + 4
 nu0 = (is-1)*9+5 ; nuf = nu0 + 4
 g0Pristine(mu0:muf,nu0:nuf) = g0Pristine(mu0:muf,nu0:nuf) + &
                                 lambda_Bi_d*LS(14:18,5:9)

 mu0 = (is-1)*9+2 ; muf = mu0 + 2
 nu0 = (is-1)*9+2 + dimH_p3; nuf = nu0 + 2
 g0Pristine(mu0:muf,nu0:nuf) = g0Pristine(mu0:muf,nu0:nuf) + &
                                 lambda_Bi_p*LS_p(1:3,4:6)
 mu0 = (is-1)*9+5 ; muf = mu0 + 4
 nu0 = (is-1)*9+5 + dimH_p3; nuf = nu0 + 4
 g0Pristine(mu0:muf,nu0:nuf) = g0Pristine(mu0:muf,nu0:nuf) + &
                                 lambda_Bi_d*LS(5:9,14:18)
 END DO

 g0i_Pristine = wc*mum_l - g0Pristine
 g0Pristine = g0i_Pristine
 CALL invers(g0Pristine,2*dimH_p3)

 g0AA = g0i_A - MATMUL(Hcll,MATMUL(g0Pristine,Hllc))
 g0BB = g0i_B - MATMUL(Hcrl,MATMUL(g0Pristine,Hrlc))

 ! saving for the next step 
 g0i_A = g0AA
 g0i_B = g0BB

 CALL invers(g0AA,dimG) 
 CALL invers(g0BB,dimG) 
 ! We have just connected flakes A and B to pristine
 ! flakes to the left and to the right, respectively.

 ! Connect flake A to the first flake of the spacer

 gAA = g0i_A - MATMUL(Hcrl,MATMUL(g0Pristine,Hrlc))
 CALL invers(gAA,dimG)

 gNN = g0i_Pristine - MATMUL(Hrlc,MATMUL(g0AA,Hcrl))
 CALL invers(gNN,2*dimH_p3)

 gAN = MATMUL(g0AA,MATMUL(Hcrl,gNN))
 gNA = MATMUL(g0Pristine,MATMUL(Hrlc,gAA))

 g0AA = gAA
 g0Nm1Nm1 = gNN
 g0ANm1 = gAN
 g0Nm1A = gNA

 ! Connect to subsequent flakes

 DO ispacer=1,Nspacer-1
 gNN = g0i_Pristine - MATMUL(Hp10,MATMUL(g0Nm1Nm1,Hp01))
 CALL invers(gNN,2*dimH_p3)
 gNA = MATMUL(gNN,MATMUL(Hp10,g0Nm1A))
 gAN = MATMUL(g0ANm1,MATMUL(Hp01,gNN))
 gAA = g0AA + MATMUL(g0ANm1,MATMUL(Hp01,gNA))
 g0AA = gAA
 g0Nm1Nm1 = gNN
 g0ANm1 = gAN
 g0Nm1A = gNA
 END DO

 ! Connect left subsystem (including adatom) to right
 ! subsystem (one flake with adatom plus one prisitne flake)
 gBB = g0i_B - MATMUL(Hcll,MATMUL(g0Nm1Nm1,Hllc))
 CALL invers(gBB,dimG)
 gBA = MATMUL(gBB,MATMUL(Hcll,g0Nm1A))
 gAB = MATMUL(g0ANm1,MATMUL(Hllc,gBB))
 gAA = g0AA + MATMUL(g0ANm1,MATMUL(Hllc,gBA))
 
 gout_ad(1,1,1:9,1:9) = gAA(1:9,1:9)
 gout_ad(1,1,1:9,10:18) = gAA(1:9,dimH+1:dimH+9)
 gout_ad(1,1,10:18,1:9) = gAA(dimH+1:dimH+9,1:9)
 gout_ad(1,1,10:18,10:18) = gAA(dimH+1:dimH+9,dimH+1:dimH+9)

 gout_ad(1,2,1:9,1:9) = gAB(1:9,1:9)
 gout_ad(1,2,1:9,10:18) = gAB(1:9,dimH+1:dimH+9)
 gout_ad(1,2,10:18,1:9) = gAB(dimH+1:dimH+9,1:9)
 gout_ad(1,2,10:18,10:18) = gAB(dimH+1:dimH+9,dimH+1:dimH+9)

 gout_ad(2,1,1:9,1:9) = gBA(1:9,1:9)
 gout_ad(2,1,1:9,10:18) = gBA(1:9,dimH+1:dimH+9)
 gout_ad(2,1,10:18,1:9) = gBA(dimH+1:dimH+9,1:9)
 gout_ad(2,1,10:18,10:18) = gBA(dimH+1:dimH+9,dimH+1:dimH+9)

 gout_ad(2,2,1:9,1:9) = gBB(1:9,1:9)
 gout_ad(2,2,1:9,10:18) = gBB(1:9,dimH+1:dimH+9)
 gout_ad(2,2,10:18,1:9) = gBB(dimH+1:dimH+9,1:9)
 gout_ad(2,2,10:18,10:18) = gBB(dimH+1:dimH+9,dimH+1:dimH+9)

 RETURN

END SUBROUTINE green

