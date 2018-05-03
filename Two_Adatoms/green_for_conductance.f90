SUBROUTINE green(wc,gNN,GNp1Np1,G01,G10)

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
 COMPLEX(double), DIMENSION(2*dimH_p3,2*dimH_p3), INTENT(OUT) :: G01,G10,gNN,gNp1Np1

 COMPLEX(double), DIMENSION(dimG,dimG) :: g0i_A,g0i_B,gaux,g0AA,g0BB,gAA,gBB,gAB,gBA
 COMPLEX(double), DIMENSION(2*dimH_p3,2*dimH_p3) :: g0i_Pristine,g0_l,gl_l,gl_r,gl_b, &
                                                    g0Spacer,g0Nm1Nm1,g0i_rl
 COMPLEX(double), DIMENSION(dimG,2*dimH_p3) :: g0ANm1,gAN
 COMPLEX(double), DIMENSION(2*dimH_p3,dimG) :: g0Nm1A,gNA

 INTEGER :: is,i,j,mu,mu0,muf,nu0,nuf,ineighb

 REAL(double) :: cri
 INTEGER :: dimdecim, itmax, ispacer

 ! Leads GF's

 ! Parameters for decim
 cri = 1d-5
 itmax = 20
 dimdecim = 2*dimH_p3

 g0i_Pristine = CMPLX(0.d0,0.d0,double)
 g0i_Pristine(1:dimH_p3,1:dimH_p3) = H00_pristine
 g0i_Pristine(dimH_p3+1:2*dimH_p3,dimH_p3+1:2*dimH_p3) = H00_pristine

 DO is=1,N_Bi
 mu0 = (is-1)*9+2; muf = mu0 + 2
 nu0 = (is-1)*9+2; nuf = nu0 + 2
 g0i_Pristine(mu0:muf,nu0:nuf) = g0i_Pristine(mu0:muf,nu0:nuf) + &
                                 lambda_Bi_p*LS_p(1:3,1:3)
 mu0 = (is-1)*9+5; muf = mu0 + 4
 nu0 = (is-1)*9+5; nuf = nu0 + 4
 g0i_Pristine(mu0:muf,nu0:nuf) = g0i_Pristine(mu0:muf,nu0:nuf) + &
                                 lambda_Bi_d*LS(5:9,5:9)

 mu0 = (is-1)*9+2 + dimH_p3; muf = mu0 + 2
 nu0 = (is-1)*9+2 + dimH_p3; nuf = nu0 + 2
 g0i_Pristine(mu0:muf,nu0:nuf) = g0i_Pristine(mu0:muf,nu0:nuf) + &
                                 lambda_Bi_p*LS_p(4:6,4:6)
 mu0 = (is-1)*9+5 + dimH_p3; muf = mu0 + 4
 nu0 = (is-1)*9+5 + dimH_p3; nuf = nu0 + 4
 g0i_Pristine(mu0:muf,nu0:nuf) = g0i_Pristine(mu0:muf,nu0:nuf) + &
                                 lambda_Bi_d*LS(14:18,14:18)

 mu0 = (is-1)*9+2 + dimH_p3; muf = mu0 + 2
 nu0 = (is-1)*9+2 ; nuf = nu0 + 2
 g0i_Pristine(mu0:muf,nu0:nuf) = g0i_Pristine(mu0:muf,nu0:nuf) + &
                                 lambda_Bi_p*LS_p(4:6,1:3)
 mu0 = (is-1)*9+5 + dimH_p3; muf = mu0 + 4
 nu0 = (is-1)*9+5 ; nuf = nu0 + 4
 g0i_Pristine(mu0:muf,nu0:nuf) = g0i_Pristine(mu0:muf,nu0:nuf) + &
                                 lambda_Bi_d*LS(14:18,5:9)

 mu0 = (is-1)*9+2 ; muf = mu0 + 2
 nu0 = (is-1)*9+2 + dimH_p3; nuf = nu0 + 2
 g0i_Pristine(mu0:muf,nu0:nuf) = g0i_Pristine(mu0:muf,nu0:nuf) + &
                                 lambda_Bi_p*LS_p(1:3,4:6)
 mu0 = (is-1)*9+5 ; muf = mu0 + 4
 nu0 = (is-1)*9+5 + dimH_p3; nuf = nu0 + 4
 g0i_Pristine(mu0:muf,nu0:nuf) = g0i_Pristine(mu0:muf,nu0:nuf) + &
                                 lambda_Bi_d*LS(5:9,14:18)
 END DO

 g0i_Pristine = wc*mum_l - g0i_Pristine
 CALL decim(g0i_Pristine,Hp01,Hp10,gl_l,gl_r,gl_b,cri,dimdecim,itmax)

 g0i_rl = gl_l 

 CALL invers(gl_l,dimdecim)
 CALL invers(gl_r,dimdecim)

 ! Flakes A amd B containing the adatoms

 g0i_A = CMPLX(0d0,0d0,double)
 g0i_B = CMPLX(0d0,0d0,double)

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

 ! Connect flake A to left lead. 
 ! The nomenclature used within the decim subroutine may
 ! be confusing: what I call the left lead is the right
 ! surface of the humongous slab generated by the decimation
 ! procedure. Accordingly, the right lead is the left surface
 ! of said slab.

 g0AA = g0i_A - MATMUL(Hcll,MATMUL(gl_r,Hllc))

 ! saving for the next step (connecting the two flakes)
 g0i_A = g0AA

 CALL invers(g0AA,dimG) 

 ! We have just connected flake A to the left
 ! leads. Now we connect the spacer 
 ! flakes to flake A, one by one.

 g0Spacer = g0i_Pristine
 CALL invers(g0Spacer,2*dimH_p3)

 ! Connect flake A to the first flake of the spacer and
 ! re-calculate the GF of flake A, the GF of the spacer 
 ! flake just added and the GF connecting both.
 gAA = g0i_A - MATMUL(Hcrl,MATMUL(g0Spacer,Hrlc))
 CALL invers(gAA,dimG)

 gNN = g0i_Pristine - MATMUL(Hrlc,MATMUL(g0AA,Hcrl))
 CALL invers(gNN,2*dimH_p3)

 gAN = MATMUL(g0AA,MATMUL(Hcrl,gNN))
 gNA = MATMUL(g0Spacer,MATMUL(Hrlc,gAA))

 g0AA = gAA
 g0Nm1Nm1 = gNN
 g0ANm1 = gAN
 g0Nm1A = gNA

 ! Now connect to subsequent flakes

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

 ! Connect flake B to the last flake of the spacer.
 gBB = g0i_B - MATMUL(Hcll,MATMUL(g0Nm1Nm1,Hllc))
 CALL invers(gBB,dimG)

 ! Now connect a pristine flake to flake B, then
 ! connect this pristine flake to the right lead. 
 ! By doing this the calculation of the conductance 
 ! is completely equivalent to that in the pristine
 ! system.

 ! Connect pristine flake
 g0iNN = g0i_Pristine - MATMUL(Hrlc,MATMUL(gBB,Hcrl))
 g0Nm1Nm1 = g0iNN
 CALL invers(g0Nm1Nm1,2*dimH_p3)

 ! Connect to right lead

 gNN = g0iNN - MATMUL(Hp01,MATMUL(gl_l,Hp10))
 CALL invers(gNN,2*dimH_p3)

 GNp1Np1 = g0i_rl - MATMUL(Hp10,MATMUL(g0Nm1Nm1,Hp01))
 CALL invers(GNp1Np1,2*dimH_p3)

 G01 = MATMUL(g0Nm1Nm1,MATMUL(Hp01,GNp1Np1))
 G10 = MATMUL(gl_l,MATMUL(Hp10,gNN))

 RETURN

END SUBROUTINE green

