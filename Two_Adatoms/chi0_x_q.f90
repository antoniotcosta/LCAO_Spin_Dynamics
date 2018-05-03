PROGRAM chi_single_TM_adatom_on_TI

 USE lattice
 USE constants_u
 USE ef_and_eta
 USE splitting
 USE sum_k
 USE lambda_SOC
 USE the_LS_matrix
 USE the_hamiltonian
 USE MPI
 USE MPI_pars


 ! Only the d part of the susceptibility matrix is
 ! relevant (U_ss and U_pp are negligible)
 COMPLEX(double), DIMENSION(4*25,25) :: chi0,chi
 COMPLEX(double), DIMENSION(4*25,4*25) :: mum,Lambda,ummLambda
 COMPLEX(double), DIMENSION(5,5) :: chipm,chipm0

 REAL(double) :: UeV,BZeeman
 REAL(double) :: m_i_A

 INTEGER :: i,j,iw,nw,is
 REAL(double) :: wmin,wmax,dw,e

 REAL(double) :: theta,phi

 REAL(double) :: ImTrchi
 COMPLEX(double) :: y
 COMPLEX(double), DIMENSION(4) :: yl ! sum rule
 COMPLEX(double), DIMENSION(4) :: sum_chi_RPA, sum_chi0

 INTEGER :: ierr,ifile
 INTEGER :: mu,nu,l,m,n

 CHARACTER(LEN=255) :: nome0
 CHARACTER(len=255), DIMENSION(4) :: nome

 REAL(double) :: qy
 REAL(double), DIMENSION(2) :: qv

  ! MPI Initialization
  CALL MPI_init(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

  pi = 4.d0*ATAN(1.d0)
  sq2 = SQRT(2.d0)
  sq3 = SQRT(3.d0)

  mum= CMPLX(0.d0,0.d0,double)
  DO mu=1,100
  mum(mu,mu) = CMPLX(1.d0,0.d0,double)
  END DO

  ef = 0d0

  OPEN(unit=99,file='chi0_entrada.in',status='old')

  ! ncp
  READ(99,*)ncp
  READ(99,*)eta
  READ(99,*)UeV 
  U = UeV
  READ(99,*)BZeeman
  READ(99,*)e0d,m_i_A
  READ(99,*)theta,phi
  READ(99,*)lambda_Bi_p,lambda_Bi_d
  READ(99,*)lambda_ad
  READ(99,*)ef
  READ(99,*)wmin,wmax,nw
  CLOSE(99)
 
  CALL hamiltonian
  CALL LS_matrix_p(theta,phi)
  CALL LS_matrix(theta,phi,LS(1,:,:))

  hw0 = BZeeman*5.7883817555d-5 ! \mu_B in eV/T

  hdel_A = 0.5d0*U*m_i_A

  IF(myid==0)THEN
  WRITE(*,*)"Input parameters:"
  WRITE(*,"('U = ',f22.12)")U
  WRITE(*,"('B = ',f22.12)")BZeeman
  WRITE(*,"('hw0 = ',f22.12)")hw0
  WRITE(*,"('e0d = ',f22.12,' m_i_A = ',f22.12)")e0D,m_i_A
  WRITE(*,"('hdel = ',f22.12)")hdel_A
  WRITE(*,"('lambda_Bi_p = ',f22.12)")lambda_Bi_p
  WRITE(*,"('lambda_Bi_d = ',f22.12)")lambda_Bi_d
  WRITE(*,"('lambda_ad = ',f22.12)")lambda_ad
  END IF

  WRITE(nome(1),"('chi_x_q.dat')")
  WRITE(nome(2),"('chi_up_x_q.dat')")
  WRITE(nome(3),"('chi_dn_x_q.dat')")
  WRITE(nome(4),"('chi_mm_x_q.dat')")


  dw = (wmax-wmin)/DBLE(nw-1)
  e = 0d0

  DO iw = 1,nw

  qv = [ wmin+dw*DBLE(iw-1), 0d0 ]
   CALL wd_Int_chi0(e,qv,chi0)
   CALL wd_Int_Lambda(e,qv,Lambda)

   ! the stuff past this point only concerns process 0,
   ! that will collect all the pieces from other processes
   ! and assemble them together to calculate chi_RPA
   IF (myid == 0) THEN

   ! Sum rule
   y = CMPLX(0d0,0d0,double)
   l=1
   DO mu=1,5
   DO nu=1,5
   m = (l-1)*25  + 5*(mu-1) + mu 
   n = 5*(nu-1) + nu
   y = y + chi0(m,n) 
   END DO
   END DO
   WRITE(*,"('Sum rule: (',D22.12,',',D22.12,')')")REAL(y),AIMAG(y)
   WRITE(888,*)e,REAL(y/(1.d0+U*y)),AIMAG(y/(1.d0+U*y))

   ! Lambdas
   yl = zero
   DO l=1,4
   DO mu=1,5
   DO nu=1,5
   m = (l-1)*25  + 5*(mu-1) + mu 
   n = 5*(nu-1) + nu
   yl(l) = yl(l) + lambda(m,n) 
   END DO
   END DO
   WRITE(*,"('Sum lambda: (',D22.12,',',D22.12,')')")REAL(yl(l)),AIMAG(yl(l))
   END DO

   ummLambda = mum - U*Lambda

   CALL invers(ummLambda,100)

   chi = MATMUL(ummLambda,chi0)
    
   sum_chi_RPA = zero
   sum_chi0 = zero

   ! Only the d part of the susceptibility 
   ! U_ss = U_pp = 0
   DO l=1,4
   DO mu=1,5
   DO nu=1,5
    m = (l-1)*25  + 5*(mu-1) + mu 
    n = 5*(nu-1) + nu
    sum_chi_RPA(l) = sum_chi_RPA(l) + chi(m,n) 
    sum_chi0(l) = sum_chi0(l) + chi0(m,n)
    IF (l==1) THEN
    chipm(mu,nu) = chi(m,n) 
    chipm0(mu,nu) = chi0(m,n) 
    END IF
  END DO
  END DO
  END DO

  OPEN(unit=789,file='Im_chi_mu_nu_x_q.dat',status='unknown',position='append')
  DO mu=1,5 
   WRITE(789,"(6(d17.8,1x))")qv(1),AIMAG(chipm(mu,:))
  END DO
  CLOSE(789)

  OPEN(unit=788,file='Im_chi0_mu_nu_x_q.dat',status='unknown',position='append')
  DO mu=1,5 
   WRITE(788,"(6(d17.8,1x))")qv(1),AIMAG(chipm0(mu,:))
  END DO
  CLOSE(788)

  OPEN(unit=789,file='Re_chi0_mu_nu_x_q.dat',status='unknown',position='append')
  DO mu=1,5 
   WRITE(789,"(6(d17.8,1x))")qv(1),REAL(chipm0(mu,:))
  END DO
  CLOSE(789)

  OPEN(unit=889,file='Re_chi_mu_nu_x_q.dat',status='unknown',position='append')
  DO mu=1,5 
   WRITE(889,"(5(d17.8,1x))")qv(1),REAL(chipm(mu,:))
  END DO
  CLOSE(889)

   
  OPEN(unit=887,file='Tr_chi_x_q.dat',status='unknown',position='append')
  ImTrchi = 0.d0
  DO mu=1,5
   ImTrchi = ImTrchi - AIMAG(chipm(mu,mu))/pi
  END DO
  WRITE(887,*)e,ImTrchi
  CLOSE(887)

  DO l=1,4
  OPEN(unit=ifile,file=nome(l),status="unknown",position="append")
  WRITE(ifile,*)qv(1),REAL(sum_chi_RPA(l)),aimag(sum_chi_RPA(l))
  CLOSE(ifile)

  WRITE(nome0,"('chi0_',I0,'_x_q.dat')")l
  OPEN(unit=ifile+1,file=nome0,status="unknown",position="append")
  WRITE(ifile+1,*)qv(1),REAL(sum_chi0(l)),aimag(sum_chi0(l))
  CLOSE(ifile+1)
  END DO

 END IF
 
 END DO   

 CALL MPI_Finalize(IERR)
 IF (IERR.ne.0) then
 WRITE(*,*)"stinky MPI finilize..."
 END IF

 STOP

END PROGRAM chi_single_TM_adatom_on_TI
