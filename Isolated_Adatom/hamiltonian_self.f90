SUBROUTINE hamiltonian_central

 USE lattice
 USE constants_u
 USE the_hamiltonian

 COMPLEX(double) :: hh

 INTEGER :: i,j,mu,nu
 REAL(double) :: Rx,Ry,Rz

 INTEGER :: numbh,ineighb,dimensao
 INTEGER :: dimH_foo

! IF (myid==0) WRITE(*,*)"dimH",dimh

 OPEN(unit=1111,file="Bismuthene.ham",status="old")

 READ(1111,*)numbh,dimH_foo
 numbneighb = numbh-1

 READ(1111,*)Rx,Ry,Rz
 DO mu=1,dimh
 DO nu=1,dimh
 READ(1111,*)hh
 H00_central(mu,nu) = CMPLX(REAL(hh),0.d0,double)
 IF(AIMAG(hh)>1d-7) WRITE(*,*)"Warning: complex hamiltonian:", i, j, mu, nu
 END DO
 END DO

 DO mu=1,dimH_foo
 H00_central(mu,mu) = H00_central(mu,mu) + deltaEF_central
 END DO

 H01_central = CMPLX(0d0,0d0,double)
 READ(1111,*)Rx,Ry,Rz
 DO mu=1,dimh
 DO nu=1,dimh
 READ(1111,*)hh
 H01_central(nu,mu) = CMPLX(REAL(hh),0.d0,double)
 IF(AIMAG(hh)>1d-7) WRITE(*,*)"Warning: complex hamiltonian:", i, j, mu, nu
 END DO
 END DO
 H01_central(dimH+1:dimG,dimH+1:dimG) = H01_central(1:dimH,1:dimH)
 H10_central = TRANSPOSE(H01_central)

 CLOSE(1111)

 Hllc = CMPLX(0d0,0d0,double)
 Hllc(1:dimH_p3,1:dimH) = H01_central(Nadatoms*Norb_ad+1:dimH,1:dimH)
 Hllc(dimH_p3+1:2*dimH_p3,dimH+1:dimG) = Hllc(1:dimH_p3,1:dimH)
 Hcll = TRANSPOSE(Hllc)

 Hcrl = CMPLX(0d0,0d0,double)
 Hcrl(1:dimH,1:dimH_p3) = H01_central(1:dimH,Nadatoms*Norb_ad+1:dimH)
 Hcrl(dimH+1:dimG,dimH_p3+1:2*dimH_p3) = Hcrl(1:dimH,1:dimH_p3)
 Hrlc = TRANSPOSE(Hcrl)


 mum = CMPLX(0d0,0d0,double)
 DO mu=1,dimG
 mum(mu,mu) = CMPLX(1d0,0d0,double)
 END DO
 

 RETURN

END SUBROUTINE hamiltonian_central

 
