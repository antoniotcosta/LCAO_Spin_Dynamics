SUBROUTINE hamiltonian_central

 USE lattice
 USE constants_u
 USE the_hamiltonian

 COMPLEX(double) :: hh

 INTEGER :: i,j,mu,nu
 REAL(double) :: Rx,Ry,Rz

 INTEGER :: numbh,ineighb,dimensao
 INTEGER :: dimH_foo

 REAL(double) :: deltaEF


 ! quick and dirty: 
 ! deltaEF = -6.109698882593924E-003
 deltaEF = -7.269101012374521D-002 ! From the original DFT calculation w/o SOC

! IF (myid==0) WRITE(*,*)"dimH",dimh

 OPEN(unit=1111,file="Bismuthene.ham",status="old")

 READ(1111,*)numbh,dimH_foo
 numbneighb = numbh-1

 WRITE(*,*)dimH,dimH_foo

 READ(1111,*)Rx,Ry,Rz
 DO mu=1,dimh
 DO nu=1,dimh
 READ(1111,*)hh
 H00_central(mu,nu) = CMPLX(REAL(hh),0.d0,double)
 IF(AIMAG(hh)>1d-7) WRITE(*,*)"Warning: complex hamiltonian:", i, j, mu, nu
 END DO
 END DO

 DO mu=1,dimH_foo
 H00_central(mu,mu) = H00_central(mu,mu) + deltaEF
 END DO

 H01_central = CMPLX(0d0,0d0,double)
 READ(1111,*)Rx,Ry,Rz
 DO mu=1,dimh
 DO nu=1,dimh
 READ(1111,*)hh
 H01_central(mu,nu) = CMPLX(REAL(hh),0.d0,double)
 IF(AIMAG(hh)>1d-7) WRITE(*,*)"Warning: complex hamiltonian:", i, j, mu, nu
 END DO
 END DO
 H01_central(dimH+1:dimG,dimH+1:dimG) = H01_central(1:dimH,1:dimH)
 H10_central = TRANSPOSE(H01_central)

 CLOSE(1111)

 mum = CMPLX(0d0,0d0,double)
 DO mu=1,dimG
 mum(mu,mu) = CMPLX(1d0,0d0,double)
 END DO
 

 RETURN

END SUBROUTINE hamiltonian_central

 
