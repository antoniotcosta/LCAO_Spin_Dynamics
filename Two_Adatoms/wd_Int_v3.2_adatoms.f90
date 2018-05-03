! In this version process 0 also has calculations to do.
! This is necessary because of SLURM's restriction on the
! number of tasks per node.
! In this version (as in the previous ones) it is assumed
! that the number of integration points is larger than the 
! number of MPI processes.

SUBROUTINE wd_Int(e,chi_hf,Lambda)

 USE lattice
 USE constants_u
 USE splitting
 USE ef_and_eta
 USE MPI_pars
 USE mpi
 USE lattice

 REAL(double), INTENT(IN) :: e
 COMPLEX(double), DIMENSION(4*5*5*Nadtotal,4*5*5*Nadtotal), INTENT(OUT) :: chi_hf,Lambda

 INTEGER, PARAMETER ::  n1=64,n2=64
 INTEGER :: i,ix,itask,send_tag
 REAL(double) :: efme,ebrk,csi,y,y2

 REAL(double), DIMENSION(n1) :: x1,p1
 REAL(double), DIMENSION(n2) :: x2,p2

 COMPLEX(double), DIMENSION(4*5*5*Nadtotal,4*5*5*Nadtotal) :: Fint,my_Fint,Lint,my_Lint

 INTEGER :: IERR,NCOUNT=4*4*5*5*5*5*Nadtotal*Nadtotal
 INTEGER, DIMENSION(MPI_STATUS_SIZE) :: stat
 
 REAL(double) :: xsgn

 LOGICAL :: isfirstpass


 CALL gauleg(0.d0,1.d0,x1,p1,n1)

 chi_hf = CMPLX(0.d0,0.d0,double)
 Lambda = CMPLX(0.d0,0.d0,double)

 IF (myid+1<=n1) THEN
  ix = myid+1
  y  = 1.d0 - x1(ix)
  y2 = y*y/( 1.d0 + eta )
  csi = ( x1(ix) + eta )/y
  WRITE(500+myid,"('Process ',I0,'  starting task ',I0,' csi = ',f21.11)")myid,ix,csi
  CALL dchi_HF_1(e,csi,my_Fint)
  my_Fint = my_Fint*p1(ix)/y2
 ELSE
  ix = myid - n1 + 1
  y  = 1.d0 - x1(ix)
  y2 = y*y/( 1.d0 + eta )
  csi = ( x1(ix) + eta )/y
  WRITE(500+myid,"('Process ',I0,'  starting task ',I0,' csi = ',f21.11)")myid,ix,csi
  CALL dLambda_1(e,csi,my_Lint)
  my_Lint = my_Lint*p1(ix)/y2
 END IF

 IF (myid==0) THEN

 ! The result of process 0's first integration point
 chi_hf = my_Fint

 DO i=1,n1-1
 CALL MPI_RECV(Fint,ncount,MPI_DOUBLE_COMPLEX,i, &
               MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)

 chi_hf = chi_hf + Fint
 WRITE(199,*)"from",stat(MPI_SOURCE),stat(MPI_TAG)
 END DO ! i=1,n1-1


 Lambda = CMPLX(0d0,0d0,double)

 DO i=n1,2*n1-1
 CALL MPI_RECV(Lint,ncount,MPI_DOUBLE_COMPLEX,i, &
               MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)

 Lambda = Lambda + Lint
 WRITE(199,*)"from",stat(MPI_SOURCE),stat(MPI_TAG)
 END DO ! i=n1,2*n1-1



 ELSE ! myid /= 0

 IF (myid+1<=n1) THEN
 CALL MPI_SEND(my_Fint,ncount,MPI_DOUBLE_COMPLEX,0,stat(MPI_TAG), &
               MPI_COMM_WORLD,ierr)
 ELSE
 CALL MPI_SEND(my_Lint,ncount,MPI_DOUBLE_COMPLEX,0,stat(MPI_TAG), &
               MPI_COMM_WORLD,ierr)

 END IF ! (myid+1<=n1)

 END IF ! myid /= 0


 IF ( ABS(e)>1.d-10 ) THEN

 IF(e>0d0) THEN
 CALL gauleg(ef-e,ef,x2,p2,n2)
 xsgn = 1d0
 ELSE
 CALL gauleg(ef,ef-e,x2,p2,n2)
 xsgn = -1d0
 END IF


 IF (myid+1<=n2) THEN
  ix = myid+1
  csi = x2(ix)
  WRITE(500+myid,"('Process ',I0,'  starting task ',I0,' csi = ',f21.11)")myid,ix,csi
  CALL dchi_HF_2(e,csi,my_Fint)
  my_Fint = my_Fint*p2(ix)
 ELSE
  ix = myid - n2 + 1
  csi = x2(ix)
  WRITE(500+myid,"('Process ',I0,'  starting task ',I0,' csi = ',f21.11)")myid,ix,csi
  CALL dLambda_2(e,csi,my_Lint)
  my_Lint = my_Lint*p2(ix)
 END IF

 IF (myid==0) THEN

 ! The result of process 0's first integration point
 chi_hf = xsgn*my_Fint

 DO i=1,n2-1
 CALL MPI_RECV(Fint,ncount,MPI_DOUBLE_COMPLEX,i, &
               MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)

 chi_hf = chi_hf + xsgn*Fint
 WRITE(199,*)"from",stat(MPI_SOURCE),stat(MPI_TAG)
 END DO ! i=1,n2-1


 Lambda = CMPLX(0d0,0d0,double)

 DO i=n2,2*n2-1
 CALL MPI_RECV(Lint,ncount,MPI_DOUBLE_COMPLEX,i, &
               MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)

 Lambda = Lambda + xsgn*Lint
 WRITE(199,*)"from",stat(MPI_SOURCE),stat(MPI_TAG)
 END DO ! i=n1,2*n1-1



 ELSE ! myid /= 0

 IF (myid+1<=n1) THEN
 CALL MPI_SEND(my_Fint,ncount,MPI_DOUBLE_COMPLEX,0,stat(MPI_TAG), &
               MPI_COMM_WORLD,ierr)
 ELSE
 CALL MPI_SEND(my_Lint,ncount,MPI_DOUBLE_COMPLEX,0,stat(MPI_TAG), &
               MPI_COMM_WORLD,ierr)

 END IF ! (myid+1<=n1)

 END IF ! myid /= 0


 END IF ! ( ABS(e)>1.d-10 )

 RETURN

END SUBROUTINE wd_Int
