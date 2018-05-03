PROGRAM Bismuthene_conductance

 USE lattice
 USE constants_u
 USE ef_and_eta
 USE lambda_SOC
 USE the_hamiltonian
 USE the_LS_matrix
 USE the_LS_p_matrix
 USE splitting

 REAL(double) :: UeV,BZeeman

 INTEGER :: i,j,nw
 REAL(double) :: wmin,wmax,dw

 REAL(double) :: theta,phi

 INTEGER ::  ileft,iright
 INTEGER :: is,js,il,ir,mu,nu,mu0,muf,nu0,nuf,ief

 REAL(double) :: Gamma0
 REAL(double), DIMENSION(N_uc-6) :: Gamma0_per_atom
 COMPLEX(double), DIMENSION(dimG,dimG) :: G00
 COMPLEX(double), DIMENSION(2*dimH_p3,2*dimH_p3) :: G11
 COMPLEX(double), DIMENSION(dimG,2*dimH_p3) :: G01,H01
 COMPLEX(double), DIMENSION(2*dimH_p3,dimG) :: G10,H10
 COMPLEX(double) :: wc

 REAL(double), DIMENSION(3,Nadtotal) :: m_i

 REAL(double) :: rho

  pi = 4.d0*ATAN(1.d0)
  sq2 = SQRT(2.d0)
  sq3 = SQRT(3.d0)

  OPEN(unit=99,file='conductance_entrada.in',status='old')
  READ(99,*)eta
  READ(99,*)UeV 
  READ(99,*)BZeeman
  READ(99,*)lambda_Bi_p,lambda_Bi_d
  READ(99,*)lambda_ad
  READ(99,*)ef,deltaEF_central,deltaEF_leads
  READ(99,*)wmin,wmax,nw
  CLOSE(99)

  OPEN(unit=99,file='self_results.in',status='old')
  DO i=1,Nadtotal
  READ(99,*)e0d(i)
  END DO
  DO i=1,Nadtotal
  DO j=1,3
  READ(99,*)m_i(j,i)
  END DO
  END DO
  CLOSE(99)

  theta = 0d0
  phi = 0d0
  hw0 = BZeeman*5.7883817555d-5 ! \mu_B in eV/T

  U = UeV
  DO i=1,Nadtotal
  hdel_A(i,1,1) = -0.5d0*U*CMPLX( m_i(3,i), 0d0,double)
  hdel_A(i,2,2) = -0.5d0*U*CMPLX(-m_i(3,i), 0d0,double)
  hdel_A(i,1,2) = -0.5d0*U*CMPLX( m_i(1,i),-m_i(2,i),double)
  hdel_A(i,2,1) = -0.5d0*U*CMPLX( m_i(1,i), m_i(2,i),double)
  END DO

  CALL hamiltonian_central
  CALL hamiltonian_leads
  CALL LS_matrix_p(theta,phi)
  CALL LS_matrix(theta,phi,LS)

  H01 = Hcrl
  H10 = Hrlc

  OPEN(unit=333,file="gamma0.dat",status="unknown",position="append")
  OPEN(unit=334,file="gamma0_per_atom.dat",status="unknown",position="append")
  OPEN(unit=666,file="DOS_central.dat",status="unknown",position="append")

  dw = (wmax-wmin)/DBLE(nw)
  DO ief = 1,nw
  ef = wmin + DBLE(ief)*dw
  wc = CMPLX(ef,eta,double)
  CALL green(wc,G00,G11,G01,G10) 
  WRITE(*,*)"sai sum"
  CALL conductance_kubo(G00,G11,G01,G10,H01,H10,Gamma0,Gamma0_per_atom)

  rho = 0d0
  DO mu=1,dimG
  rho = rho - AIMAG(G00(mu,mu))/pi
  END DO 
  WRITE(666,*)ef,rho

  WRITE(333,*)ef,Gamma0
  WRITE(334,"(110(f19.10,1x))")ef,Gamma0_per_atom

  END DO

  CLOSE(333)
  CLOSE(334)
  CLOSE(666)

  STOP

END PROGRAM Bismuthene_conductance



