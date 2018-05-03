! In this version process 0 also has calculations to do.
! This is necessary because of SLURM's restriction on the
! number of tasks per node.
! In this version (as in the previous ones) it is assumed
! that the number of integration points is larger than the 
! number of MPI processes.

SUBROUTINE wd_Int_chi0(e,chi_hf)

 USE lattice
 USE constants_u
 USE splitting
 USE ef_and_eta
 USE MPI_pars
 USE mpi
 USE lattice

 REAL(double), INTENT(IN) :: e
 COMPLEX(double), DIMENSION(4*5*5*Nadtotal,4*5*5*Nadtotal), INTENT(OUT) :: chi_hf

 INTEGER, PARAMETER ::  n1=64,n2=64
 INTEGER :: i,ix,itask,send_tag
 REAL(double) :: efme,ebrk,csi,y,y2

 REAL(double), DIMENSION(n1) :: x1,p1
 REAL(double), DIMENSION(n2) :: x2,p2

 COMPLEX(double), DIMENSION(4*5*5*Nadtotal,4*5*5*Nadtotal) :: Fint,my_Fint

 INTEGER :: IERR,NCOUNT=4*4*5*5*5*5*Nadtotal*Nadtotal
 INTEGER, DIMENSION(MPI_STATUS_SIZE) :: stat
 
 REAL(double) :: xsgn

 LOGICAL :: isfirstpass


 CALL gauleg(0.d0,1.d0,x1,p1,n1)

 chi_hf = CMPLX(0.d0,0.d0,double)

 IF (myid==0) THEN

 ! initialize all processes
 DO i=1,numprocs-1
 send_tag = i
 CALL MPI_SEND(i,1,MPI_INTEGER,i,send_tag,MPI_COMM_WORLD,ierr)
 END DO

 ! This is process 0's first point 
 ix = numprocs
 y  = 1.d0 - x1(ix)
 y2 = y*y/( 1.d0 + eta )
 csi = ( x1(ix) + eta )/y

 WRITE(50,"('Process 0 starting task ',I0,' csi = ',f21.11)"),ix,csi
 CALL dchi_HF_1(e,csi,my_Fint)

 WRITE(50,*)"Process 0 finished task",ix
 my_Fint = my_Fint*p1(ix)/y2

 itask = numprocs + 1

 ELSE ! myid /= 0
 CALL MPI_RECV(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
 WRITE(50+myid,*)"sou",myid,"recebi tarefa inicial"
 END IF

 IF (myid==0) THEN

 ! The result of process 0's first integration point
 chi_hf = my_Fint

 isfirstpass = .TRUE.

 DO WHILE (itask<=n1.OR.isfirstpass)

 isfirstpass = .FALSE.

 DO i=1,numprocs-1
 CALL MPI_RECV(Fint,ncount,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE, &
               MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)

 chi_hf = chi_hf + Fint
 WRITE(199,*)"from",stat(MPI_SOURCE),stat(MPI_TAG)

 IF (itask<=n1) THEN
 CALL MPI_SEND(itask,1,MPI_INTEGER,stat(MPI_SOURCE),itask, &
               MPI_COMM_WORLD,ierr)
 itask = itask + 1
 ELSE
 CALL MPI_SEND(0,1,MPI_INTEGER,stat(MPI_SOURCE),0,MPI_COMM_WORLD,ierr)
 END IF

 END DO ! i=1,numprocs-1

 ! If there is any remainder integration point after all other
 ! processes got new points, process 0 will take it
 IF (i>=numprocs.AND.itask<=n1) THEN
 WRITE(50,*)"Getting remainder task",itask
 ix = itask
 y  = 1.d0 - x1(ix)
 y2 = y*y/( 1.d0 + eta )
 csi = ( x1(ix) + eta )/y

 CALL dchi_HF_1(e,csi,my_Fint)
 WRITE(50,*)"Finished remainder task"

 chi_hf = chi_hf + my_Fint*p1(ix)/y2

 itask = itask + 1
 END IF

 END DO ! WHILE (itask<=n1)

 ELSE ! myid /= 0

 DO WHILE (stat(MPI_TAG)/=0)

 y  = 1.d0 - x1(ix)
 y2 = y*y/( 1.d0 + eta )
 csi = ( x1(ix) + eta )/y
 WRITE(50+myid,"('Starting calculation of task ',I0,' csi = ',f21.11)"),ix,csi
 CALL dchi_HF_1(e,csi,Fint)

 WRITE(50+myid,*)"Finished calculation of task",ix
 Fint = Fint*p1(ix)/y2

 CALL MPI_SEND(Fint,ncount,MPI_DOUBLE_COMPLEX,0,stat(MPI_TAG), &
               MPI_COMM_WORLD,ierr)

 CALL MPI_RECV(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD, &
               stat,ierr)

 WRITE(50+myid,*)"sou",myid,"recebi",stat(MPI_TAG),ix

 END DO

 END IF


 IF ( ABS(e)>1.d-10 ) THEN

 IF(e>0d0) THEN
 CALL gauleg(ef-e,ef,x2,p2,n2)
 xsgn = 1d0
 ELSE
 CALL gauleg(ef,ef-e,x2,p2,n2)
 xsgn = -1d0
 END IF

 IF (myid==0) THEN

 ! initialize all processes
 DO i=1,numprocs-1
 send_tag = i
 CALL MPI_SEND(i,1,MPI_INTEGER,i,send_tag,MPI_COMM_WORLD,ierr)
 END DO

 ! This is process 0's first point 
 ix = numprocs
 csi = x2(ix)
 CALL dchi_HF_2(e,csi,my_Fint)

 my_Fint = my_Fint*p2(ix)

 itask = numprocs + 1

 ELSE
 
 CALL MPI_RECV(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
 ! WRITE(50+myid,*)"sou",myid,"recebi tarefa inicial"

 END IF ! (myid==0)

 IF (myid==0) THEN

 ! The result of process 0's first integration point
 chi_hf = chi_hf + xsgn*my_Fint

 isfirstpass = .TRUE.

 DO WHILE (itask<=n2.OR.isfirstpass)

 isfirstpass = .FALSE.

 DO i=1,numprocs-1
 CALL MPI_RECV(Fint,ncount,MPI_DOUBLE_COMPLEX, &
               MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)

 chi_hf = chi_hf + xsgn*Fint

 ! WRITE(97,*)"from",stat(MPI_SOURCE),stat(MPI_TAG)

 IF (itask<=n2) THEN
 CALL MPI_SEND(itask,1,MPI_INTEGER,stat(MPI_SOURCE), &
               itask,MPI_COMM_WORLD,ierr)

 itask = itask + 1
 ELSE
 ! Finalize tasks
 CALL MPI_SEND(0,1,MPI_INTEGER,stat(MPI_SOURCE),0,MPI_COMM_WORLD,ierr)
 END IF

 END DO ! i=1,numprocs-1

 ! If there is any remainder integration point after all other
 ! processes got new points, process 0 will take it
 IF (i>=numprocs.AND.itask<=n2) THEN
 ix = itask
 csi = x2(ix)
 CALL dchi_HF_2(e,csi,my_Fint)
 chi_hf = chi_hf + xsgn*my_Fint*p2(ix)
 itask = itask + 1
 END IF

 END DO ! WHILE (itask<=n2)

 ELSE

 DO WHILE (stat(MPI_TAG)/=0)

 csi = x2(ix)

 CALL dchi_HF_2(e,csi,Fint)

 Fint = Fint*p2(ix)

 CALL MPI_SEND(Fint,ncount,MPI_DOUBLE_COMPLEX,0,stat(MPI_TAG),MPI_COMM_WORLD,ierr)

 CALL MPI_RECV(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)

 ! WRITE(50+myid,*)"sou",myid,"recebi",stat(MPI_TAG),ix

 END DO

 END IF

 END IF

 RETURN

END SUBROUTINE wd_Int_chi0 





SUBROUTINE wd_Int_Lambda(e,Lambda)

 USE lattice
 USE constants_u
 USE splitting
 USE ef_and_eta
 USE MPI_pars
 USE mpi
 USE mpi_constants
 USE lattice

 REAL(double), INTENT(IN) :: e
 COMPLEX(double), DIMENSION(4*5*5*Nadtotal,4*5*5*Nadtotal), INTENT(OUT) :: Lambda

 INTEGER, PARAMETER ::  n1=64,n2=64
 INTEGER :: i,ix,itask,send_tag
 REAL(double) :: efme,ebrk,csi,y,y2

 REAL(double), DIMENSION(n1) :: x1,p1
 REAL(double), DIMENSION(n2) :: x2,p2

 COMPLEX(double), DIMENSION(4*5*5*Nadtotal,4*5*5*Nadtotal) :: Fint,my_Fint

 INTEGER :: IERR,NCOUNT=4*4*5*5*5*5*Nadtotal*Nadtotal
 INTEGER, DIMENSION(MPI_STATUS_SIZE) :: stat

 REAL(double) :: xsgn

 LOGICAL :: isfirstpass

 CALL gauleg(0.d0,1.d0,x1,p1,n1)

 Lambda = CMPLX(0.d0,0.d0,double)

 IF (myid==0) THEN

 ! initialize all processes
 DO i=1,numprocs-1
 send_tag = i
 CALL MPI_SEND(i,1,MPI_INTEGER,i,send_tag,MPI_COMM_WORLD,ierr)
 END DO

 ! This is process 0's first point 
 ix = numprocs
 y  = 1.d0 - x1(ix)
 y2 = y*y/( 1.d0 + eta )
 csi = ( x1(ix) + eta )/y

 CALL dLambda_1(e,csi,my_Fint)

 my_Fint = my_Fint*p1(ix)/y2
 itask = numprocs + 1

 ELSE ! myid /= 0
 CALL MPI_RECV(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
 ! WRITE(50+myid,*)"sou",myid,"recebi tarefa inicial"
 END IF

 IF (myid==0) THEN

 ! The result of process 0's first integration point
 Lambda = my_Fint

 isfirstpass = .TRUE.

 DO WHILE (itask<=n1.OR.isfirstpass)

 isfirstpass = .FALSE.

 DO i=1,numprocs-1
 CALL MPI_RECV(Fint,ncount,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE, &
               MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)

 Lambda = Lambda + Fint
 ! WRITE(199,*)"from",stat(MPI_SOURCE),stat(MPI_TAG)

 IF (itask<=n1) THEN
 CALL MPI_SEND(itask,1,MPI_INTEGER,stat(MPI_SOURCE),itask,MPI_COMM_WORLD,ierr)

 itask = itask + 1
 ELSE
 CALL MPI_SEND(0,1,MPI_INTEGER,stat(MPI_SOURCE),0,MPI_COMM_WORLD,ierr)
 END IF

 END DO ! i=1,numprocs-1

 ! If there is any remainder integration point after all other
 ! processes got new points, process 0 will take it
 IF (i>=numprocs.AND.itask<=n1) THEN
 WRITE(50,*)"Getting remainder task",itask
 ix = itask
 y  = 1.d0 - x1(ix)
 y2 = y*y/( 1.d0 + eta )
 csi = ( x1(ix) + eta )/y

 CALL dLambda_1(e,csi,my_Fint)
 WRITE(50,*)"Finished remainder task"

 Lambda = Lambda + my_Fint*p1(ix)/y2

 itask = itask + 1
 END IF

 END DO ! WHILE (itask<=n1)

 ELSE ! myid /= 0

 DO WHILE (stat(MPI_TAG)/=0)

 y  = 1.d0 - x1(ix)
 y2 = y*y/( 1.d0 + eta )
 csi = ( x1(ix) + eta )/y

 CALL dLambda_1(e,csi,Fint)

 Fint = Fint*p1(ix)/y2

 CALL MPI_SEND(Fint,ncount,MPI_DOUBLE_COMPLEX,0,stat(MPI_TAG),MPI_COMM_WORLD,ierr)

 CALL MPI_RECV(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)

 ! WRITE(50+myid,*)"sou",myid,"recebi",stat(MPI_TAG),ix

 END DO

 END IF


 IF ( ABS(e)>1.d-10 ) THEN

 IF(e>0d0) THEN
 CALL gauleg(ef-e,ef,x2,p2,n2)
 xsgn = 1d0
 ELSE
 CALL gauleg(ef,ef-e,x2,p2,n2)
 xsgn = -1d0
 END IF 

 IF (myid==0) THEN

 ! initialize all processes
 DO i=1,numprocs-1
 send_tag = i
 CALL MPI_SEND(i,1,MPI_INTEGER,i,send_tag,MPI_COMM_WORLD,ierr)
 END DO

 ! This is process 0's first point 
 ix = numprocs
 csi = x2(ix)
 CALL dLambda_2(e,csi,my_Fint)
 my_Fint = my_Fint*p2(ix)

 itask = numprocs + 1

 ELSE
 
 CALL MPI_RECV(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
 ! WRITE(50+myid,*)"sou",myid,"recebi tarefa inicial"

 END IF ! (myid==0)

 IF (myid==0) THEN
   
 ! The result of process 0's first integration point
 Lambda = Lambda + xsgn*my_Fint

 isfirstpass = .TRUE.

 DO WHILE (itask<=n2.OR.isfirstpass)

 isfirstpass = .FALSE.

 DO i=1,numprocs-1
 CALL MPI_RECV(Fint,ncount,MPI_DOUBLE_COMPLEX, &
               MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)

 Lambda = Lambda + xsgn*Fint

 ! WRITE(97,*)"from",stat(MPI_SOURCE),stat(MPI_TAG)

 IF (itask<=n2) THEN
 CALL MPI_SEND(itask,1,MPI_INTEGER,stat(MPI_SOURCE),itask,MPI_COMM_WORLD,ierr)

 itask = itask + 1
 ELSE
 ! Finalize tasks
 CALL MPI_SEND(0,1,MPI_INTEGER,stat(MPI_SOURCE),0,MPI_COMM_WORLD,ierr)
 END IF

 END DO ! i=1,numprocs-1

 ! If there is any remainder integration point after all other
 ! processes got new points, process 0 will take it
 IF (i>=numprocs.AND.itask<=n2) THEN
 ix = itask
 csi = x2(ix)
 CALL dLambda_2(e,csi,my_Fint)
 Lambda = Lambda + xsgn*my_Fint*p2(ix)
 itask = itask + 1
 END IF

 END DO ! WHILE (itask<=n2)

 ELSE

 DO WHILE (stat(MPI_TAG)/=0)

 csi = x2(ix)

 CALL dLambda_2(e,csi,Fint)

 Fint = Fint*p2(ix)

 CALL MPI_SEND(Fint,ncount,MPI_DOUBLE_COMPLEX,0,stat(MPI_TAG),MPI_COMM_WORLD,ierr)

 CALL MPI_RECV(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)

 ! WRITE(50+myid,*)"sou",myid,"recebi",stat(MPI_TAG),ix

 END DO

 END IF

 END IF ! ABS(e) > 1d-10

 RETURN

END SUBROUTINE wd_Int_Lambda 

