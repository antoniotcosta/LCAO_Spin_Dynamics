PROGRAM DOS_graphene_with_adatom

 USE lattice
 USE constants_u
 USE ef_and_eta
 USE lambda_SOC
 USE splitting
 USE the_LS_p_matrix
 USE the_LS_matrix
 USE the_hamiltonian, ONLY: deltaEF_central,deltaEF_leads


 INTEGER :: i,j,iw,nw,is,iTM,mu,nu
 REAL(double) :: wmin,wmax,dw

 COMPLEX(double) :: wc

 REAL(double), DIMENSION(Nadtotal) :: rho_Fe_up,rho_Fe_dn
 COMPLEX(double) :: GAB_up,GBA_up,GAB_dn,GBA_dn

 COMPLEX(double), DIMENSION(Nadtotal,Nadtotal,18,18) :: G_ad

 REAL(double) :: theta0,phi0, UeV, BZeeman

 REAL(double), DIMENSION(3,Nadtotal) :: m_i

 CHARACTER(LEN=255) :: fmtstrng,nome

  pi = 4.d0*ATAN(1.d0)
  sq2 = SQRT(2.d0)
  sq3 = SQRT(3.d0)

  OPEN(unit=99,file='DOS_entrada.in',status='old')

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
 
  theta0 = 0d0
  phi0 = 0d0
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
  CALL LS_matrix_p(theta0,phi0)
  CALL LS_matrix(theta0,phi0,LS)


  OPEN(unit=101,file="DOS_Fe_1.dat",status='unknown',position='append')
  OPEN(unit=102,file="DOS_Fe_2.dat",status='unknown',position='append')

  dw = (wmax-wmin)/DBLE(nw-1)
  DO iw=1,nw

   wc = CMPLX(wmin + DBLE(iw-1)*dw, eta)
   CALL green(wc,G_ad)

   rho_Fe_up = 0d0
   rho_Fe_dn = 0d0
   DO is=1,Nadtotal
   DO mu=1,9
   rho_Fe_up(is) = rho_Fe_up(is)-AIMAG(G_ad(is,is,mu,mu))/pi
   rho_Fe_dn(is) = rho_Fe_dn(is)-AIMAG(G_ad(is,is,mu+9,mu+9))/pi
   END DO
   END DO
   WRITE(101,*)REAL(wc),rho_Fe_up(1),rho_Fe_dn(1)
   WRITE(102,*)REAL(wc),rho_Fe_up(2),rho_Fe_dn(2)

   gAB_up = CMPLX(0d0,0d0,double)
   gAB_dn = CMPLX(0d0,0d0,double)
   gBA_up = CMPLX(0d0,0d0,double)
   gBA_dn = CMPLX(0d0,0d0,double)
   DO mu=1,9
   DO nu=1,9
   gAB_up = gAB_up + G_ad(1,2,mu,nu)
   gAB_dn = gAB_dn + G_ad(1,2,mu+9,nu+9)
   gBA_up = gBA_up + G_ad(2,1,mu,nu)
   gBA_dn = gBA_dn + G_ad(2,1,mu+9,nu+9)
   END DO
   END DO
   WRITE(201,*)REAL(wc),REAL(gAB_up), AIMAG(gAB_up)
   WRITE(202,*)REAL(wc),REAL(gAB_dn), AIMAG(gAB_dn)
   WRITE(203,*)REAL(wc),REAL(gBA_up), AIMAG(gBA_up)
   WRITE(204,*)REAL(wc),REAL(gBA_dn), AIMAG(gBA_dn)

  END DO

  DO is=1,Nadtotal
  CLOSE(100+is)
  END DO

  STOP

END PROGRAM DOS_graphene_with_adatom


