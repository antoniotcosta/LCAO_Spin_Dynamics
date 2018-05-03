INCLUDE '/home/antc/lib/f90_kind.f90'
INCLUDE '/home/antc/lib/constants_u.f90'

MODULE  lattice

 USE f90_kind

 INTEGER, PARAMETER :: N_Bi = 108, Nadatoms = 1, N_uc = N_Bi+NAdatoms+6 
 INTEGER, PARAMETER :: Norb_Bi = 9, Norb_ad = 9
 INTEGER :: dimh 
 INTEGER :: dimG
 INTEGER, PARAMETER :: nnmax=8
 REAL(double), DIMENSION(3) :: b1,b2
 REAL(double), DIMENSION(N_uc,3) :: r0

 INTEGER :: numbneighb ! Read in routine hamiltonian

 REAL(double), DIMENSION(nnmax,3) :: dij

END MODULE lattice


MODULE the_hamiltonian

 USE lattice

 REAL(double), DIMENSION(:,:), ALLOCATABLE, SAVE :: H00
 REAL(double), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: Hij

END MODULE the_hamiltonian


MODULE ef_and_eta

 USE f90_kind

 REAL(double) :: ef,eta

END MODULE ef_and_eta



MODULE sum_k

 USE f90_kind

 INTEGER :: ncp

END MODULE sum_k


MODULE the_LS_matrix

 USE lattice

 COMPLEX(double), DIMENSION(18,18) :: LS

END MODULE the_LS_matrix

MODULE the_LS_p_matrix

 USE lattice

 COMPLEX(double), DIMENSION(6,6) :: LS_p

END MODULE the_LS_p_matrix


MODULE lambda_SOC

 USE lattice

 REAL(double) :: lambda_Bi_p,lambda_Bi_d,lambda_ad

END MODULE lambda_SOC


MODULE splitting

 USE lattice

 REAL(double), DIMENSION(Nadatoms) :: hdel_A,e0d
 REAL(double) :: hw0,U

END MODULE splitting



!--------------------------End of MODULES section--------------------------



PROGRAM DOS_graphene_with_N_adatoms

 USE lattice
 USE constants_u
 USE ef_and_eta
 USE sum_k
 USE lambda_SOC
 USE the_LS_matrix
 USE splitting
 USE H_of_k

 REAL(double) :: tol,xtol
 REAL(double) :: UeV,BZeeman

 INTEGER :: i,j,iq,iw,is,iTM,mu,mu0,muf
 REAL(double) :: qmin,qmax,dq

 CHARACTER(LEN=255) :: fmtstrng,nome

 INTEGER :: IERR

 CHARACTER(len=1) :: option

 INTEGER :: n,lda, ifail
 INTEGER :: lwork 
 CHARACTER(LEN=1) :: job,uplo
 REAL(double), DIMENSION(:), ALLOCATABLE :: ek
 REAL(double), DIMENSION(:), ALLOCATABLE :: rwork 
 REAL(double), DIMENSION(:), ALLOCATABLE :: work 

 REAL(double), DIMENSION(2) :: qv 
 REAL(double) :: kx

 REAL(double) :: psi_up,psi_dn
 REAL(double) :: theta0,phi0

 REAL(double), DIMENSION(Nadatoms) :: m_i_A

 COMPLEX(double), DIMENSION(:,:), ALLOCATABLE :: gout


  pi = 4.d0*ATAN(1.d0)
  sq2 = SQRT(2.d0)
  sq3 = SQRT(3.d0)

  OPEN(unit=99,file='bands_entrada.in',status='old')

  ! ncp
  READ(99,*)ncp
  ! Coulomb interaction (1 eV=1/13.6 Ry) 
  READ(99,*)lambda_Bi_p,lambda_Bi_d
  READ(99,*)lambda_ad
  READ(99,*)U
  READ(99,*)e0d
  READ(99,*)m_i_A
  READ(99,*)qmin,qmax

  CLOSE(99)

  theta0 = 0d0
  phi0 = 0d0

  CALL LS_matrix_p(theta0,phi0)
  CALL LS_matrix(theta0,phi0,LS)
  CALL hamiltonian

  hdel_A = 0.5d0*U*m_i_A

  dimG = 2*dimH
  lwork = 6*dimG
  ALLOCATE(gout(dimG,dimG))
  ALLOCATE(ek(dimG))
  ALLOCATE(rwork(10*dimG))
  ALLOCATE(work(lwork))

  WRITE(*,*)"eigenvaLue or eigenVector?"
  READ(*,*)option

  IF (option=="L") THEN

  dq = (qmax-qmin)/DBLE(ncp-1)

  DO iq=0,ncp-1

  qv = [qmin+dq*DBLE(iq),0.d0]

  CALL Hk(qv,gout)

  n = dimG
  job = 'V'  
  uplo = 'U'
  lda = n
  ifail = 0
  CALL f02haf(job,uplo,n,gout,lda,ek,rwork,work,lwork,ifail)
 
  DO i=1,dimG
  IF (ABS(ek(i))<2d0) WRITE(800+i,"(9(f15.6,1x))")qv(1),ek(i)
  END DO

  END DO

  DO i=1,dimG
  CLOSE(800+i)
  END DO

  ELSE ! IF(OPTION=="V") THEN

  WRITE(*,*)"Wavevector?"
  READ(*,*)kx
  qv = [kx,0.d0]
  CALL Hk(qv,gout)

  n = dimG
  job = 'V'  
  uplo = 'U'
  lda = n
  ifail = 0
  CALL f02haf(job,uplo,n,gout,lda,ek,rwork,work,lwork,ifail)

  DO iw=1,dimG
  WRITE(nome,"('eigen_',I0,'.dat')")iw
  OPEN(unit=1800,file=nome,status='unknown') !,position='rewind')
  DO i=1,N_uc
  psi_up = SUM(gout((i-1)*8+1:(i-1)*8+4,iw)*CONJG(gout((i-1)*8+1:(i-1)*8+4,iw)))
  psi_dn = SUM(gout((i-1)*8+5:(i-1)*8+8,iw)*CONJG(gout((i-1)*8+5:(i-1)*8+8,iw)))
  WRITE(1800,"(I0,2(f18.9,1x))")i,psi_up,psi_dn 
  END DO
  
	  CLOSE(1800)
  END DO

  END IF 

  STOP

END PROGRAM DOS_graphene_with_N_adatoms
