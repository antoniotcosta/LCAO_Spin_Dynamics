PROGRAM Self_Bismuthene_with_adatoms

 USE lattice
 USE constants_u
 USE ef_and_eta
 USE param_fcn
 USE splitting
 USE lambda_SOC
 USE the_LS_matrix
 USE the_LS_p_matrix
 USE the_hamiltonian, ONLY : deltaEF_leads,deltaEF_central
 USE MPI_pars
 USE MPI

 REAL(double) :: tol,xtol
 REAL(double) :: UeV,BZeeman

 INTEGER :: i,j,iw,nw,is,iTM,mu
 REAL(double) :: wmin,wmax,dw
 REAL(double) :: lambda_C_val

 REAL(double), DIMENSION(3,Nadtotal) :: m_i_A

 CHARACTER(LEN=255) :: fmtstrng,nome

 INTEGER :: IERR

 INTEGER :: ifail,NX
 INTEGER, PARAMETER :: N=4*Nadtotal,lwa=N*(3*N+13)/2    !
 REAL(double), DIMENSION(4*Nadtotal) :: X,fvec
 REAL(double), DIMENSION(lwa) :: wa             !
 
 REAL(double) :: theta0,phi0

 EXTERNAL fcn


  ! MPI Initialization
  CALL MPI_init(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

  pi = 4.d0*ATAN(1.d0)
  sq2 = SQRT(2.d0)
  sq3 = SQRT(3.d0)

  OPEN(unit=99,file='self_entrada.in',status='old')

  ! ncp
  READ(99,*)eta
  ! Coulomb interaction (1 eV=1/13.6 Ry) 
  READ(99,*)UeV ! = 1.d0
  U = UeV
  READ(99,*)BZeeman
  READ(99,*)lambda_Bi_p,Lambda_Bi_d
  READ(99,*)lambda_ad
  READ(99,*)xnofel
  READ(99,*)ef,deltaEF_central,deltaEF_leads
  READ(99,*)Nspacer

  CLOSE(99)

  OPEN(unit=99,file='e0d_miA.in',status='old')
  DO i=1,Nadtotal
  READ(99,*)e0d(i),m_i_A(:,i)
  END DO
  CLOSE(99)

  hw0 = BZeeman*5.7883817555d-5 ! \mu_B in eV/T

  theta0 = 0d0
  phi0 = 0d0

  CALL LS_matrix_p(theta0,phi0)
  CALL LS_matrix(theta0,phi0,LS)

  CALL hamiltonian_leads
  CALL hamiltonian_central

  tol  = 1.d-9
  xtol = 1.d-9
  NX = 4*Nadtotal

  X(1:Nadtotal) = e0d
  DO i=1,Nadtotal
  X(Nadtotal+1+(i-1)*3:Nadtotal+1+(i-1)*3+2) = m_i_A(:,i)
  END DO
 
  IF(myid==0) THEN
   WRITE(*,*)"Input parameters:"
   DO i=1,Nadtotal
   WRITE(*,"('e0d(',I0,') = ',f15.5)")i,e0d(i)
   WRITE(*,"('m0(',I0,') = ',3(f15.5,1x))")i,m_i_A(:,i)
   END DO
   WRITE(*,"('U = ',f15.5)")U
   WRITE(*,"('EF = ',f15.5)")ef
   WRITE(*,"('Nspacer = ',I0)")Nspacer
  END IF


  ifail = -1 ! noisy exit
  CALL c05nbf(fcn,NX,X,fvec,xtol,wa,lwa,ifail)
  !CALL fcn(NX,X,fvec,ifail)
  !ifail = 0
  IF (ifail/=0) THEN
  WRITE(*,*)"c05nbf failed; terminating..."
  CALL MPI_Finalize(IERR)
  IF (IERR.ne.0) then
     WRITE(*,*)"stinky MPI finilize..."
  END IF
  STOP
  END IF

  IF (myid==0) THEN
  WRITE(*,*)
  WRITE(*,*)"Finished successfully."
  WRITE(*,*)

  DO i=1,4*Nadtotal
  WRITE(222,*)X(i)
  END DO

  END IF


  CALL MPI_Finalize(IERR)
  IF (IERR.ne.0) then
     WRITE(*,*)"stinky MPI finilize..."
  END IF

  STOP

END PROGRAM Self_Bismuthene_with_adatoms


SUBROUTINE fcn(NX,X,fvec,iflag)

 USE lattice
 USE constants_u
 USE param_fcn
 USE ef_and_eta
 USE splitting
 USE the_LS_matrix
 USE MPI
 USE MPI_pars

 INTEGER, INTENT(IN) :: NX,iflag
 REAL(double), DIMENSION(NX), INTENT(IN) :: X
 REAL(double), DIMENSION(NX), INTENT(OUT) :: fvec

 INTEGER :: n_o_nproc,resto,acum_n
 INTEGER, DIMENSION(0:99) :: npppr
 INTEGER :: ninte,i,mu
 INTEGER, PARAMETER :: nptint = 64
 REAL(double) :: y,weoy2
 REAL(double) :: xe(nptint),we(nptint)
 REAL(double) :: my_xe(0:99,nptint),my_we(0:99,nptint)
 REAL(double), DIMENSION(Nadtotal,Norb_ad) :: n_orb_u,n_orb_d,dn_orb_u,dn_orb_d,dSx,dSy
 REAL(double), DIMENSION(Nadtotal,Norb_ad) :: my_n_orb_u,my_n_orb_d
 REAL(double), DIMENSION(Nadtotal,Norb_ad) :: their_n_orb_u,their_n_orb_d
 REAL(double), DIMENSION(Nadtotal) :: my_Sx,my_Sy,their_Sx,their_Sy,Sx,Sy
 COMPLEX(double) :: ec
 COMPLEX(double), DIMENSION(Nadtotal,Nadtotal,2*Norb_ad,2*Norb_ad) :: gout_ad

 REAL(double), DIMENSION(3,Nadtotal) :: m_i,m_f
 REAL(double), DIMENSION(Nadtotal) ::  n_t
 
 INTEGER, DIMENSION(MPI_STATUS_SIZE) :: stat
 INTEGER :: ierr

 INTEGER :: pksz = Nadtotal*Norb_ad

 INTEGER :: il,ix,is

 

  e0d = X(1:Nadtotal)
  DO i=1,Nadtotal
  m_i(:,i) = X(Nadtotal+1+(i-1)*3:Nadtotal+1+(i-1)*3+2)
  END DO
  DO i=1,Nadtotal
  hdel_A(i,1,1) = -0.5d0*U*CMPLX( m_i(3,i), 0d0,double)
  hdel_A(i,2,2) = -0.5d0*U*CMPLX(-m_i(3,i), 0d0,double)
  hdel_A(i,1,2) = -0.5d0*U*CMPLX( m_i(1,i),-m_i(2,i),double)
  hdel_A(i,2,1) = -0.5d0*U*CMPLX( m_i(1,i), m_i(2,i),double)
  END DO

  IF (myid==0) THEN
  WRITE(*,*)"------Input to fcn------"
  DO i=1,Nadtotal
  WRITE(*,"('e0d(',I0,') = ',(f15.5))")i,e0d(i)
  WRITE(*,"('m_i(',I0,') = ',3(f15.5,1x))")i,m_i(:,i)
  END DO
  DO i=1,Nadtotal
  WRITE(*,"('------hdel(',I0,')------')")i
  WRITE(*,"(2(2(f15.5,1x),5x))")REAL(hdel_A(i,1,1)),AIMAG(hdel_A(i,1,1)),REAL(hdel_A(i,1,2)),AIMAG(hdel_A(i,1,2))
  WRITE(*,"(2(2(f15.5,1x),5x))")REAL(hdel_A(i,2,1)),AIMAG(hdel_A(i,2,1)),REAL(hdel_A(i,2,2)),AIMAG(hdel_A(i,2,2))
  END DO
  END IF

  ! Integration in energy

  ninte = nptint                               
  CALL gauleg(0.d0,1.d0,xe,we,ninte)

  n_o_nproc = ninte/numprocs
  ! numero de pontos por processo
  resto = MOD(ninte,numprocs)
  IF (resto>0) npppr(0:resto-1) = n_o_nproc + 1
  npppr(resto:numprocs-1) = n_o_nproc

  my_xe(0,1:npppr(0)) = xe(1:npppr(0))
  my_we(0,1:npppr(0)) = we(1:npppr(0))
 
  DO i=1,numprocs-1
  acum_n = SUM(npppr(0:i-1))
  my_xe(i,1:npppr(i)) = xe( acum_n+1 : acum_n+npppr(i) )
  my_we(i,1:npppr(i)) = we( acum_n+1 : acum_n+npppr(i) )
  END DO

  my_n_orb_u = 0.d0
  my_n_orb_d = 0.d0
  my_Sx = 0.d0
  my_Sy = 0.d0
  

  DO i=1,npppr(myid)
   y  = 1.d0 - my_xe(myid,i)
   ec = CMPLX(ef,(my_xe(myid,i)+eta)/y)
   weoy2 = (1.d0+eta)*my_we(myid,i)/(y*y)

   CALL green(ec,gout_ad)

   ! Adatoms occupancies
   DO is =1,Nadtotal
   FORALL (mu=1:Norb_ad)
   dn_orb_u(is,mu) = REAL(gout_ad(is,is,mu,mu))
   dn_orb_d(is,mu) = REAL(gout_ad(is,is,mu+Norb_ad,mu+Norb_ad))
   dSx(is,mu) =  REAL(gout_ad(is,is,mu+Norb_ad,mu)) +  &
                             REAL(gout_ad(is,is,mu,mu+Norb_ad))
   dSy(is,mu) =  AIMAG(gout_ad(is,is,mu+Norb_ad,mu)) -  &
                            AIMAG(gout_ad(is,is,mu,mu+Norb_ad)) 
   END FORALL
   END DO

   my_n_orb_u = my_n_orb_u + weoy2*dn_orb_u
   my_n_orb_d = my_n_orb_d + weoy2*dn_orb_d
   DO is=1,Nadtotal
   my_Sx(is) = my_Sx(is) + SUM(dSx(is,Norb_ad-4:Norb_ad))*weoy2
   my_Sy(is) = my_Sy(is) + SUM(dSy(is,Norb_ad-4:Norb_ad))*weoy2
   END DO

  END DO

  ! Collect and add up
  IF (myid == 0) THEN
   n_orb_u = my_n_orb_u
   n_orb_d = my_n_orb_d
   Sx = my_Sx
   Sy = my_Sy

   DO i=1,numprocs-1
    CALL MPI_Recv(their_n_orb_u,pksz,MPI_double_precision, &
                  i,1000+i,MPI_comm_world,stat,ierr) 
    CALL MPI_Recv(their_n_orb_d,pksz,MPI_double_precision, &
                  i,2000+i,MPI_comm_world,stat,ierr) 
    CALL MPI_Recv(their_Sx,Nadtotal,MPI_double_precision,      &
                      i,3000+i,MPI_comm_world,stat,ierr) 
    CALL MPI_Recv(their_Sy,Nadtotal,MPI_double_precision,      &
                      i,4000+i,MPI_comm_world,stat,ierr) 

    n_orb_u = n_orb_u + their_n_orb_u
    n_orb_d = n_orb_d + their_n_orb_d
    Sx = Sx + their_Sx
    Sy = Sy + their_Sy

         
   END DO
  ELSE
   CALL MPI_send(my_n_orb_u,pksz,MPI_double_precision,0,  &
                 1000+myid,MPI_comm_world,ierr)
   CALL MPI_send(my_n_orb_d,pksz,MPI_double_precision,0,  &
                 2000+myid,MPI_comm_world,ierr)
   CALL MPI_send(my_Sx,Nadtotal,MPI_double_precision,0,       &
                   3000+myid,MPI_comm_world,ierr)
   CALL MPI_send(my_Sy,Nadtotal,MPI_double_precision,0,       &
                   4000+myid,MPI_comm_world,ierr)

  END IF


  CALL MPI_bcast(n_orb_u,pksz,MPI_double_precision,0,MPI_COMM_WORLD,ierr)
  CALL MPI_bcast(n_orb_d,pksz,MPI_double_precision,0,MPI_COMM_WORLD,ierr)
  CALL MPI_bcast(Sx,Nadtotal,MPI_double_precision,0,MPI_COMM_WORLD,ierr)
  CALL MPI_bcast(Sy,Nadtotal,MPI_double_precision,0,MPI_COMM_WORLD,ierr)

  n_orb_u = n_orb_u/pi + .5d0
  n_orb_d = n_orb_d/pi + .5d0
  Sx = Sx/pi
  Sy = Sy/pi

  DO is=1,Nadtotal
  n_t(is) = SUM( n_orb_u(is,:) + n_orb_d(is,:) )
  m_f(1,is) = Sx(is)
  m_f(2,is) = Sy(is)
  m_f(3,is) = SUM(n_orb_u(is,Norb_ad-4:Norb_ad) - n_orb_d(is,Norb_ad-4:Norb_ad))
  END DO

  DO is=1,Nadtotal
  fvec(is) =  xnofel(is) - SUM(n_orb_u(is,5:9)+n_orb_d(is,5:9))
  fvec(Nadtotal+(is-1)*3+1) = m_f(1,is) - m_i(1,is)
  fvec(Nadtotal+(is-1)*3+2) =  m_f(2,is) - m_i(2,is)
  fvec(Nadtotal+(is-1)*3+3) =  m_f(3,is) - m_i(3,is)
  END DO

  IF (myid==0) THEN
   WRITE(*,"('--------output----------')")
   DO is=1,Nadtotal
   WRITE(*,"('n_T(',I0,') = ',e18.8)")is,n_t(is)
   WRITE(*,"('n_d(',I0,') = ',e18.8)")is,SUM(n_orb_u(is,5:9)+n_orb_d(is,5:9))
   WRITE(*,"('delta_n_d(',I0,') = ',e18.8)")is,fvec(is)
   WRITE(*,"('m_f(',I0,') = ',3(e15.5,1x))")is,m_f(:,is)
   WRITE(*,"('delta_m(',I0,') = ',3(e15.5,1x))")is,fvec(Nadtotal+(is-1)*3+1: &
                                                        Nadtotal+(is-1)*3+3)

   WRITE(*,*)
   END DO
  END IF

  RETURN


END SUBROUTINE fcn








