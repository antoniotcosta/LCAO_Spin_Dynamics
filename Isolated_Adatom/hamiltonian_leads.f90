SUBROUTINE hamiltonian_leads

 USE lattice
 USE constants_u
 USE the_hamiltonian
! USE MPI_pars

 COMPLEX(double) :: hh

 INTEGER :: i,j,mu,nu,mu0,muf,nu0,nuf
 REAL(double) :: Rx,Ry,Rz

 INTEGER :: numbh,ineighb,dimensao
 INTEGER :: dimH_foo

! IF (myid==0) WRITE(*,*)"dimH",dimh

 OPEN(unit=1111,file="Bismuthene_pristine_3x1.ham",status="old")

 READ(1111,*)numbh,dimH_foo
 numbneighb = numbh-1

 READ(1111,*)Rx,Ry,Rz
 DO mu=1,dimh_p3
 DO nu=1,dimh_p3
 READ(1111,*)hh
 H00_pristine(mu,nu) = CMPLX(REAL(hh),0.d0,double)
 IF(AIMAG(hh)>1d-7) WRITE(*,*)"Warning: complex hamiltonian:", i, j, mu, nu
 END DO
 END DO

 DO ineighb = 1,numbneighb

 READ(1111,*)dij(ineighb,:)
 DO mu=1,dimh_p3
 DO nu=1,dimh_p3
 READ(1111,*)hh
 H01_pristine(nu,mu) = CMPLX(REAL(hh),0.d0,double)
 IF(AIMAG(hh)>1d-7) WRITE(*,*)"Warning: complex hamiltonian:", i, j, mu, nu
 END DO
 END DO

 END DO

 ! Aligning EF to zero
 DO mu=1,dimh_p3
 H00_pristine(mu,mu) = H00_pristine(mu,mu) + deltaEF_leads ! .169336d0
 END DO

 CLOSE(1111)

 tl = CMPLX(0d0,0d0,double)
 tld = CMPLX(0d0,0d0,double)
 tl(1:dimH_p3,1:dimH_p3) = H01_pristine
 tl(dimH_p3+1:2*dimH_p3,dimH_p3+1:2*dimH_p3) = tl(1:dimH_p3,1:dimH_p3)
 tld(1:dimH_p3,1:dimH_p3) = TRANSPOSE(H01_pristine)
 tld(dimH_p3+1:2*dimH_p3,dimH_p3+1:2*dimH_p3) = tld(1:dimH_p3,1:dimH_p3)


 mum_l = CMPLX(0d0,0d0,double)
 DO i=1,2*dimH_p3
 mum_l(i,i) = CMPLX(1d0,0d0,double) 
 END DO

 RETURN

END SUBROUTINE hamiltonian_leads
