SUBROUTINE hamiltonian_central

 USE lattice
 USE constants_u
 USE the_hamiltonian

 COMPLEX(double) :: hh

 INTEGER :: i,j,mu,nu
 REAL(double) :: Rx,Ry,Rz

 INTEGER :: numbh,ineighb,dimensao
 INTEGER :: dimH_foo

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

 HAB = CMPLX(0d0,0d0,double)
 READ(1111,*)Rx,Ry,Rz
 DO mu=1,dimh
 DO nu=1,dimh
 READ(1111,*)hh
 HAB(nu,mu) = CMPLX(REAL(hh),0.d0,double)
 IF(AIMAG(hh)>1d-7) WRITE(*,*)"Warning: complex hamiltonian:", i, j, mu, nu
 END DO
 END DO
 HAB(dimH+1:dimG,dimH+1:dimG) = HAB(1:dimH,1:dimH)
 HBA = TRANSPOSE(HAB)

 CLOSE(1111)

 mum = CMPLX(0d0,0d0,double)
 DO i=1,dimG
 mum(i,i) = CMPLX(1d0,0d0,double) 
 END DO


 RETURN

END SUBROUTINE hamiltonian_central

 
