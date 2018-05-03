SUBROUTINE hamiltonian_central

 USE lattice
 USE constants_u
 USE the_hamiltonian
! USE MPI_pars

 COMPLEX(double) :: hh

 INTEGER :: i,j,mu,nu
 REAL(double) :: Rx,Ry,Rz

 INTEGER :: numbh,ineighb,dimensao
 INTEGER :: dimH_foo

! IF (myid==0) WRITE(*,*)"dimH",dimh

 OPEN(unit=1111,file="Bismuthene.ham",status="old")

 READ(1111,*)numbh,dimH
 numbneighb = numbh-1

 ALLOCATE(H00(dimH,dimH))
 ALLOCATE(Hij(numbneighb,dimH,dimH))

 READ(1111,*)Rx,Ry,Rz
 DO mu=1,dimh
 DO nu=1,dimh
 READ(1111,*)hh
 H00(mu,nu) = CMPLX(REAL(hh),0.d0,double)
 IF(AIMAG(hh)>1d-7) WRITE(*,*)"Warning: complex hamiltonian:", i, j, mu, nu
 END DO
 END DO

 DO ineighb = 1,numbneighb

 READ(1111,*)dij(ineighb,:)
 DO mu=1,dimh
 DO nu=1,dimh
 READ(1111,*)hh
 Hij(ineighb,mu,nu) = CMPLX(REAL(hh),0.d0,double)
 IF(AIMAG(hh)>1d-7) WRITE(*,*)"Warning: complex hamiltonian:", i, j, mu, nu
 END DO
 END DO

 END DO

 CLOSE(1111)

 RETURN

END SUBROUTINE hamiltonian_central

 
