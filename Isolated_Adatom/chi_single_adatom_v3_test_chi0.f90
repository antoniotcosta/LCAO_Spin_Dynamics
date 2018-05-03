PROGRAM chi_single_TM_adatom_on_TI

 USE lattice
 USE constants_u
 USE ef_and_eta
 USE splitting
 USE lambda_SOC
 USE the_LS_matrix
 USE the_hamiltonian
 USE MPI
 USE MPI_pars


 ! Only the d part of the susceptibility matrix is
 ! relevant (U_ss and U_pp are negligible)
 COMPLEX(double), DIMENSION(4*25*Nadtotal,4*25*Nadtotal) :: chi0,chi
 COMPLEX(double), DIMENSION(4*25*Nadtotal,4*25*Nadtotal) :: mum_chi,Lambda,ummLambda

 REAL(double) :: UeV,BZeeman

 INTEGER :: i,j,iw,nw
 REAL(double) :: wmin,wmax,dw,e

 REAL(double) :: theta,phi

 REAL(double) :: ImTrchi
 COMPLEX(double) :: y
 COMPLEX(double), DIMENSION(4,4,Nadtotal,Nadtotal) :: sum_chi_RPA, sum_chi0


 REAL(double), DIMENSION(3,Nadtotal) :: m_i

 INTEGER :: ierr,ifile
 INTEGER :: mu,nu,l,ll,m,n,is,js

 CHARACTER(LEN=255) :: nome0
 CHARACTER(len=255) :: nome

  ! MPI Initialization
  CALL MPI_init(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

  pi = 4.d0*ATAN(1.d0)
  sq2 = SQRT(2.d0)
  sq3 = SQRT(3.d0)

  mum_chi= CMPLX(0.d0,0.d0,double)
  DO mu=1,100*Nadtotal
  mum_chi(mu,mu) = CMPLX(1.d0,0.d0,double)
  END DO

  OPEN(unit=99,file='chi_entrada.in',status='old')

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

  CALL hamiltonian_leads
  CALL hamiltonian_central
  CALL LS_matrix(theta,phi,LS)
  CALL LS_matrix_p(theta,phi)


  IF(myid==0)THEN
  WRITE(*,*)"Input parameters:"
  WRITE(*,"('U = ',f22.12)")U
  WRITE(*,"('B = ',f22.12)")BZeeman
  WRITE(*,"('hw0 = ',f22.12)")hw0
  DO i=1,Nadtotal
  WRITE(*,"('e0d(',I0,') = ',f22.12,' m_i(',I0,') = ',3(f22.12,1x))")i,e0D(i),i,m_i(:,i)
  END DO
  ! WRITE(*,"('hdel = ',f22.12)")hdel_A
  WRITE(*,"('lambda_Bi_p = ',f22.12)")lambda_Bi_p
  WRITE(*,"('lambda_Bi_d = ',f22.12)")lambda_Bi_d
  WRITE(*,"('lambda_ad = ',f22.12)")lambda_ad
  END IF

  IF (nw>1) THEN
  dw = (wmax-wmin)/DBLE(nw-1)
  ELSE
  dw = 0d0
  END IF

  DO iw = 1,nw

   e = wmin+dw*DBLE(iw-1)

   CALL wd_Int_chi0(e,chi0)
   ! CALL wd_Int_Lambda(e,Lambda)

   ! the stuff past this point only concerns process 0,
   ! that will collect all the pieces from other processes
   ! and assemble them together to calculate chi_RPA
   IF (myid == 0) THEN

   ummLambda = mum_chi - U*Lambda

   CALL invers(ummLambda,100*Nadtotal)

   chi = MATMUL(ummLambda,chi0)
    
   sum_chi_RPA = zero
   sum_chi0 = zero

   ! Only the d part of the susceptibility 
   ! U_ss = U_pp = 0
   DO l=1,4
   DO ll=1,4
   DO is=1,Nadtotal
   DO js=1,Nadtotal
   DO mu=1,5
   DO nu=1,5
    m = (l-1)*25*Nadtotal + (is-1)*25 + (mu-1)*5 + mu 
    n = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(nu-1) + nu
    sum_chi_RPA(l,ll,is,js) = sum_chi_RPA(l,ll,is,js) + chi(m,n) 
    sum_chi0(l,ll,is,js) = sum_chi0(l,ll,is,js) + chi0(m,n)
   END DO
   END DO
   END DO
   END DO
   END DO
   END DO

  DO l=1,4
  DO ll=1,4
  DO is=1,Nadtotal
  DO js=1,Nadtotal
  WRITE(nome,"('chi_',I0,'_',I0,'_',I0,'_',I0,'.dat')")l,ll,is,js
  OPEN(unit=999,file=nome,status="unknown",position="append")
  WRITE(999,*)e,REAL(sum_chi_RPA(l,ll,is,js)),aimag(sum_chi_RPA(l,ll,is,js))
  CLOSE(999)
  END DO
  END DO
  END DO
  END DO

  DO l=1,4
  DO ll=1,4
  DO is=1,Nadtotal
  DO js=1,Nadtotal
  WRITE(nome,"('chi0_',I0,'_',I0,'_',I0,'_',I0,'.dat')")l,ll,is,js
  OPEN(unit=999,file=nome,status="unknown",position="append")
  WRITE(999,*)e,REAL(sum_chi0(l,ll,is,js)),aimag(sum_chi0(l,ll,is,js))
  CLOSE(999)
  END DO
  END DO
  END DO
  END DO

 END IF
 
 END DO   

 CALL MPI_Finalize(IERR)
 IF (IERR.ne.0) then
 WRITE(*,*)"stinky MPI finilize..."
 END IF

 STOP

END PROGRAM chi_single_TM_adatom_on_TI
