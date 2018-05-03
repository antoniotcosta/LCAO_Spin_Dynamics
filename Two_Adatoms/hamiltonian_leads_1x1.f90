SUBROUTINE hamiltonian_leads

 USE lattice
 USE constants_u
 USE the_hamiltonian
! USE MPI_pars

 COMPLEX(double) :: hh

 COMPLEX(double), DIMENSION(dimh_p1-2,dimh_p1-2) :: Haux

 INTEGER :: i,j,mu,nu,mu0,muf,nu0,nuf
 REAL(double) :: Rx,Ry,Rz

 INTEGER :: numbh,ineighb,dimensao
 INTEGER :: dimH_foo

! IF (myid==0) WRITE(*,*)"dimH",dimh

 OPEN(unit=1111,file="Bismuthene_pristine.ham",status="old")

 READ(1111,*)numbh,dimH_foo
 numbneighb = numbh-1

 READ(1111,*)Rx,Ry,Rz
 DO mu=1,dimh_p1
 DO nu=1,dimh_p1
 READ(1111,*)hh
 H00(mu,nu) = CMPLX(REAL(hh),0.d0,double)
 IF(AIMAG(hh)>1d-7) WRITE(*,*)"Warning: complex hamiltonian:", i, j, mu, nu
 END DO
 END DO

 DO ineighb = 1,numbneighb

 READ(1111,*)dij(ineighb,:)
 DO mu=1,dimh_p1
 DO nu=1,dimh_p1
 READ(1111,*)hh
 Hij(ineighb,mu,nu) = CMPLX(REAL(hh),0.d0,double)
 IF(AIMAG(hh)>1d-7) WRITE(*,*)"Warning: complex hamiltonian:", i, j, mu, nu
 END DO
 END DO

 END DO

 CLOSE(1111)

 ! Building tripled cell
 ! push H to the end of the matrix
 ! Original positions for H atoms: 1 and 20
 dimHBi = dimH_p1 - 2

 H00_pristine = CMPLX(0d0,0d0,double)

 ! Bi: 2-19 (18 atoms) + 21-38 (18 atoms)
 H00_pristine(1:162,1:162) = H00(2:163,2:163)
 H00_pristine(163:324,163:324) = H00(165:326,165:326)
 H00_pristine(1:162,163:324) = H00(2:163,165:326)
 H00_pristine(163:324,1:162) = H00(165:326,2:163)

 H00_pristine(325:648,325:648) = H00_pristine(1:324,1:324)
 H00_pristine(649:972,649:972) = H00_pristine(1:324,1:324)


 Haux(  1:162,  1:162) = Hij(1,  2:163,  2:163)
 Haux(163:324,163:324) = Hij(1,165:326,165:326)
 Haux(  1:162,163:324) = Hij(1,  2:163,165:326)
 Haux(163:324,  1:162) = Hij(1,165:326,  2:163)

 H00_pristine(1:324,325:648) = Haux 
 H00_pristine(325:648,1:324) = TRANSPOSE(Haux)

 H00_pristine(325:648,649:972) = Haux 
 H00_pristine(649:972,325:648) = TRANSPOSE(Haux)

 Haux(  1:162,  1:162) = Hij(2,  2:163,  2:163)
 Haux(163:324,163:324) = Hij(2,165:326,165:326)
 Haux(  1:162,163:324) = Hij(2,  2:163,165:326)
 Haux(163:324,  1:162) = Hij(2,165:326,  2:163)

 H00_pristine(1:324,649:972) = Haux 
 H00_pristine(649:972,1:324) = TRANSPOSE(Haux)

 ! H atoms
 H00_pristine(973,973) = H00(1,1)
 H00_pristine(974,974) = H00(164,164)
 H00_pristine(975,975) = H00(1,1)
 H00_pristine(976,976) = H00(164,164)
 H00_pristine(977,977) = H00(1,1)
 H00_pristine(978,978) = H00(164,164)

 H00_pristine(973,  1:162) = H00(1,  2:163)
 H00_pristine(973,163:324) = H00(1,165:326)

 H00_pristine(974,  1:162) = H00(164,  2:163)
 H00_pristine(974,163:324) = H00(164,165:326)

 H00_pristine(975,325:486) = H00(1,  2:163)
 H00_pristine(975,487:648) = H00(1,165:326)

 H00_pristine(976,325:486) = H00(164,  2:163)
 H00_pristine(976,487:648) = H00(164,165:326)

 H00_pristine(977,649:810) = H00(1,  2:163)
 H00_pristine(977,811:972) = H00(1,165:326)

 H00_pristine(978,649:810) = H00(164,  2:163)
 H00_pristine(978,811:972) = H00(164,165:326)

 H00_pristine(  1:162,973) = H00(2:163,1)
 H00_pristine(163:324,973) = H00(165:326,1)

 H00_pristine(  1:162,974) = H00(  2:163,164)
 H00_pristine(163:324,974) = H00(165:326,164)

 H00_pristine(325:486,975) = H00(2:163,1)
 H00_pristine(487:648,975) = H00(1,165:326)

 H00_pristine(325:486,976) = H00(  2:163,164)
 H00_pristine(487:648,976) = H00(165:326,164) 

 H00_pristine(649:810,977) = H00(2:163,1)
 H00_pristine(811:972,977) = H00(165:326,1)

 H00_pristine(649:810,978) = H00(  2:163,164)
 H00_pristine(811:972,978) = H00(165:326,164)

 ! Inter-cell
 H01_pristine = CMPLX(0d0,0d0,double)

 Haux(  1:162,  1:162) = Hij(3,  2:163,  2:163)
 Haux(163:324,163:324) = Hij(3,165:326,165:326)
 Haux(  1:162,163:324) = Hij(3,  2:163,165:326)
 Haux(163:324,  1:162) = Hij(3,165:326,  2:163)

 DO i=1,3
 mu0 = (i-1)*dimHBi + 1; muf = mu0 + dimHBi - 1
 H01_pristine(mu0:muf,mu0:muf) = Haux
 END DO

 Haux(  1:162,  1:162) = Hij(2,  2:163,  2:163)
 Haux(163:324,163:324) = Hij(2,165:326,165:326)
 Haux(  1:162,163:324) = Hij(2,  2:163,165:326)
 Haux(163:324,  1:162) = Hij(2,165:326,  2:163)

 i=2; j=1
 mu0 = (i-1)*dimHBi + 1; muf = mu0 + dimHBi - 1
 nu0 = (j-1)*dimHBi + 1; nuf = nu0 + dimHBi - 1
 H01_pristine(mu0:muf,nu0:nuf) = Haux

 i=3; j=2
 mu0 = (i-1)*dimHBi + 1; muf = mu0 + dimHBi - 1
 nu0 = (j-1)*dimHBi + 1; nuf = nu0 + dimHBi - 1
 H01_pristine(mu0:muf,nu0:nuf) = Haux

 Haux(  1:162,  1:162) = Hij(1,  2:163,  2:163)
 Haux(163:324,163:324) = Hij(1,165:326,165:326)
 Haux(  1:162,163:324) = Hij(1,  2:163,165:326)
 Haux(163:324,  1:162) = Hij(1,165:326,  2:163)

 i=3; j=1
 mu0 = (i-1)*dimHBi + 1; muf = mu0 + dimHBi - 1
 nu0 = (j-1)*dimHBi + 1; nuf = nu0 + dimHBi - 1
 H01_pristine(mu0:muf,nu0:nuf) = Haux

 tl = CMPLX(0d0,0d0,double)
 tld = CMPLX(0d0,0d0,double)
 tl(1:dimH_p3,1:dimH_p3) = H01_pristine
 tl(dimH_p3+1:2*dimH_p3,dimH_p3+1:2*dimH_p3) = tl(1:dimH_p3,1:dimH_p3)
 tld(1:dimH_p3,1:dimH_p3) = TRANSPOSE(H01_pristine)
 tld(dimH_p3+1:2*dimH_p3,dimH_p3+1:2*dimH_p3) = tld(1:dimH_p3,1:dimH_p3)


 Hllc = CMPLX(0d0,0d0,double)
 Hcll = CMPLX(0d0,0d0,double)
 Hllc(1:2*dimH_p3,2*Nadatoms*Norb_ad+1:dimG) = tl
 Hcll(2*Nadatoms*Norb_ad+1:dimG,1:2*dimH_p3) = tld

 Hrlc = CMPLX(0d0,0d0,double)
 Hcrl = CMPLX(0d0,0d0,double)
 Hrlc(1:2*dimH_p3,2*Nadatoms*Norb_ad+1:dimG) = tld
 Hcrl(2*Nadatoms*Norb_ad+1:dimG,1:2*dimH_p3) = tl

 mum = CMPLX(0d0,0d0,double)
 DO i=1,dimG
 mum(i,i) = CMPLX(1d0,0d0,double) 
 END DO

 mum_l = CMPLX(0d0,0d0,double)
 DO i=1,2*dimH_p3
 mum_l(i,i) = CMPLX(1d0,0d0,double) 
 END DO

 RETURN

END SUBROUTINE hamiltonian_leads
