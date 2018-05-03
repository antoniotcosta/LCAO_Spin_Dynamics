INCLUDE '/home/antc/lib/f90_kind.f90'
INCLUDE '/home/antc/lib/constants_u.f90'



MODULE  lattice

 USE f90_kind

 INTEGER, PARAMETER :: N_Bi = 36, Nadatoms = 1, N_uc = N_Bi+NAdatoms+2 
 INTEGER, PARAMETER :: Norb_Bi = 9, Norb_ad = 9
 INTEGER, PARAMETER :: dimh = N_Bi*Norb_Bi + Nadatoms*Norb_ad + 2
 INTEGER, PARAMETER :: dimG = 2*dimH
 INTEGER, PARAMETER :: nnmax=8
 REAL(double), DIMENSION(3) :: b1,b2
 REAL(double), DIMENSION(N_uc,3) :: r0

 INTEGER :: numbneighb ! Read in routine hamiltonian

 REAL(double), DIMENSION(nnmax,3) :: dij

END MODULE lattice



MODULE the_hamiltonian

 USE lattice

 REAL(double), DIMENSION(dimH,dimH) :: H00
 REAL(double), DIMENSION(nnmax,dimH,dimH) :: Hij

END MODULE the_hamiltonian



module splitting

 USE lattice

 REAL(double), DIMENSION(Nadatoms) :: hdel_A,e0d
 REAL(double) :: hw0,U

end module splitting



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

 COMPLEX(double), DIMENSION(Nadatoms,18,18) :: LS

END MODULE the_LS_matrix



MODULE the_LS_p_matrix

 USE lattice

 COMPLEX(double), DIMENSION(6,6) :: LS_p

END MODULE the_LS_p_matrix

MODULE lambda_SOC

 USE lattice

 REAL(double) :: lambda_Bi_p,lambda_Bi_d,lambda_ad

END MODULE lambda_SOC



MODULE param_fcn

 USE lattice

 REAL(double), DIMENSION(Nadatoms) :: xnofel

END MODULE param_fcn
 

MODULE MPI_pars

 USE f90_kind

 INTEGER :: myid,numprocs
   
END MODULE MPI_pars



!--------------------------End of MODULES section--------------------------



PROGRAM Self_Bismuthene_with_adatoms

 USE lattice
 USE constants_u
 USE ef_and_eta
 USE param_fcn
 USE splitting
 USE sum_k
 USE lambda_SOC
 USE the_LS_p_matrix
 USE the_LS_matrix
 USE MPI_pars
 USE MPI

 REAL(double) :: tol,xtol
 REAL(double) :: UeV,BZeeman

 INTEGER :: i,j,iw,nw,is,iTM,mu
 REAL(double) :: wmin,wmax,dw
 REAL(double) :: lambda_C_val

 REAL(double), DIMENSION(Nadatoms) :: theta,phi,m_i_A

 CHARACTER(LEN=255) :: fmtstrng,nome

 INTEGER :: IERR

 INTEGER :: ifail,NX
 INTEGER, PARAMETER :: N=4*Nadatoms,lwa=N*(3*N+13)/2    !
 REAL(double), DIMENSION(4*Nadatoms) :: X,fvec
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

  OPEN(unit=99,file='find_EF_entrada.in',status='old')

  ! ncp
  READ(99,*)ncp
  READ(99,*)eta
  ! Coulomb interaction (1 eV=1/13.6 Ry) 
  READ(99,*)UeV ! = 1.d0
  U = UeV
  READ(99,*)BZeeman
  READ(99,*)lambda_Bi_p,Lambda_Bi_d
  READ(99,*)lambda_ad
  READ(99,*)ef,xnofel

  CLOSE(99)

  hw0 = 0d0
  hdel_A = 0d0
  e0d = 0d0

  CALL hamiltonian

  theta0 = 0d0
  phi0 = 0d0
  CALL LS_matrix(theta0,phi0,LS(1,:,:))
  CALL LS_matrix_p(theta0,phi0)

  CALL fcn

  !tol  = 1.d-9
  !xtol = 1.d-9
  !NX = 4*Nadatoms

  !X(1:Nadatoms) = e0d
  !X(Nadatoms+1:2*Nadatoms) = m_i_A
  !X(2*Nadatoms+1:3*Nadatoms) = theta
  !X(3*Nadatoms+1:4*Nadatoms) = phi
 
  !IF(myid==0) WRITE(*,*)X
  !IF(myid==0) WRITE(*,*)U,ef

  !ifail = -1 ! noisy exit
  !CALL c05nbf(fcn,NX,X,fvec,xtol,wa,lwa,ifail)

  !IF (ifail/=0) THEN
  !WRITE(*,*)"c05nbf failed; terminating..."
  !CALL MPI_Finalize(IERR)
  !IF (IERR.ne.0) then
  !   WRITE(*,*)"stinky MPI finilize..."
  !END IF
  !STOP
  !END IF


  CALL MPI_Finalize(IERR)
  IF (IERR.ne.0) then
     WRITE(*,*)"stinky MPI finilize..."
  END IF

  STOP

END PROGRAM Self_Bismuthene_with_adatoms


!SUBROUTINE fcn(NX,X,fvec,iflag)
SUBROUTINE fcn

 USE lattice
 USE constants_u
 USE param_fcn
 USE ef_and_eta
 USE splitting
 USE the_LS_matrix
 USE the_LS_p_matrix
 USE MPI
 USE MPI_pars

! INTEGER, INTENT(IN) :: NX,iflag
! REAL(double), DIMENSION(NX), INTENT(IN) :: X
! REAL(double), DIMENSION(NX), INTENT(OUT) :: fvec

 INTEGER :: n_o_nproc,resto,acum_n
 INTEGER, DIMENSION(0:99) :: npppr
 INTEGER :: ninte,i,mu
 INTEGER, PARAMETER :: nptint = 64
 REAL(double) :: y,weoy2
 REAL(double) :: xe(nptint),we(nptint)
 REAL(double) :: my_xe(0:99,nptint),my_we(0:99,nptint)
 REAL(double), DIMENSION(Nadatoms,Norb_ad) :: n_orb_u,n_orb_d,dn_orb_u,dn_orb_d,dSx,dSy
 REAL(double), DIMENSION(Nadatoms,Norb_ad) :: my_n_orb_u,my_n_orb_d
 REAL(double), DIMENSION(Nadatoms,Norb_ad) :: their_n_orb_u,their_n_orb_d
 REAL(double), DIMENSION(Nadatoms) :: my_Sx,my_Sy,their_Sx,their_Sy,Sx,Sy
 REAL(double), DIMENSION(N_Bi,Norb_Bi) :: n_orb_u_Bi,n_orb_d_Bi,&
                                          dn_orb_u_Bi,dn_orb_d_Bi
 REAL(double), DIMENSION(N_Bi,Norb_Bi) :: my_n_orb_u_Bi,my_n_orb_d_Bi
 REAL(double), DIMENSION(N_Bi,Norb_Bi) :: their_n_orb_u_Bi,their_n_orb_d_Bi
 COMPLEX(double) :: ec
 COMPLEX(double), DIMENSION(Nadatoms,Nadatoms,2*Norb_ad,2*Norb_ad) :: gout_ad
 COMPLEX(double), DIMENSION(N_Bi,N_Bi,2*Norb_Bi,2*Norb_Bi) :: gout_Bi

 REAL(double), DIMENSION(Nadatoms) :: m_i,theta,phi,n_t,m_f
 
 INTEGER, DIMENSION(MPI_STATUS_SIZE) :: stat
 INTEGER :: ierr

 INTEGER :: pksz = Nadatoms*Norb_ad, pksz_Bi = N_Bi*Norb_Bi

 INTEGER :: il,ix,is

 REAL(double) :: n_up_Bi,n_dn_Bi,n_T_Bi
 

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
  my_n_orb_u_Bi = 0.d0 
  my_n_orb_d_Bi = 0.d0 
  my_Sx = 0.d0
  my_Sy = 0.d0
  

  DO i=1,npppr(myid)
   y  = 1.d0 - my_xe(myid,i)
   ec = CMPLX(ef,(my_xe(myid,i)+eta)/y)
   weoy2 = (1.d0+eta)*my_we(myid,i)/(y*y)

   CALL sumC_linear(ec,gout_Bi,gout_ad)

   ! C occupancies
   DO is =1,N_Bi
   FORALL (mu=1:Norb_Bi)
   dn_orb_u_Bi(is,mu) = REAL(gout_Bi(is,is,mu,mu))
   dn_orb_d_Bi(is,mu) = REAL(gout_Bi(is,is,mu+Norb_Bi,mu+Norb_Bi))
   END FORALL
   END DO

   ! Adatoms occupancies
   DO is =1,Nadatoms
   FORALL (mu=1:Norb_ad)
   dn_orb_u(is,mu) = REAL(gout_ad(is,is,mu,mu))
   dn_orb_d(is,mu) = REAL(gout_ad(is,is,mu+Norb_ad,mu+Norb_ad))
   dSx(is,mu) = 0.5d0*(  REAL(gout_ad(is,is,mu+Norb_ad,mu)) +  &
                             REAL(gout_ad(is,is,mu,mu+Norb_ad)) ) 
   dSy(is,mu) = 0.5d0*( AIMAG(gout_ad(is,is,mu+Norb_ad,mu)) -  &
                            AIMAG(gout_ad(is,is,mu,mu+Norb_ad)) ) 
   END FORALL
   END DO

   my_n_orb_u_Bi = my_n_orb_u_Bi + weoy2*dn_orb_u_Bi
   my_n_orb_d_Bi = my_n_orb_d_Bi + weoy2*dn_orb_d_Bi
   my_n_orb_u = my_n_orb_u + weoy2*dn_orb_u
   my_n_orb_d = my_n_orb_d + weoy2*dn_orb_d
   DO is=1,Nadatoms
   my_Sx(is) = my_Sx(is) + SUM(dSx(is,Norb_ad-4:Norb_ad))*weoy2
   my_Sy(is) = my_Sy(is) + SUM(dSy(is,Norb_ad-4:Norb_ad))*weoy2
   END DO

  END DO

  ! Collect and add up
  IF (myid == 0) THEN
   n_orb_u = my_n_orb_u
   n_orb_d = my_n_orb_d
   n_orb_u_Bi = my_n_orb_d_Bi
   n_orb_d_Bi = my_n_orb_d_Bi
   Sx = my_Sx
   Sy = my_Sy

   DO i=1,numprocs-1
    CALL MPI_Recv(their_n_orb_u_Bi,pksz_Bi,MPI_double_precision, &
                  i,5000+i,MPI_comm_world,stat,ierr) 
    CALL MPI_Recv(their_n_orb_d_Bi,pksz_Bi,MPI_double_precision, &
                  i,6000+i,MPI_comm_world,stat,ierr) 
    CALL MPI_Recv(their_n_orb_u,pksz,MPI_double_precision, &
                  i,1000+i,MPI_comm_world,stat,ierr) 
    CALL MPI_Recv(their_n_orb_d,pksz,MPI_double_precision, &
                  i,2000+i,MPI_comm_world,stat,ierr) 
    CALL MPI_Recv(their_Sx,Nadatoms,MPI_double_precision,      &
                      i,3000+i,MPI_comm_world,stat,ierr) 
    CALL MPI_Recv(their_Sy,Nadatoms,MPI_double_precision,      &
                      i,4000+i,MPI_comm_world,stat,ierr) 

    n_orb_u_Bi = n_orb_u_Bi + their_n_orb_u_Bi
    n_orb_d_Bi = n_orb_d_Bi + their_n_orb_d_Bi
    n_orb_u = n_orb_u + their_n_orb_u
    n_orb_d = n_orb_d + their_n_orb_d
    Sx = Sx + their_Sx
    Sy = Sy + their_Sy

         
   END DO
  ELSE
   CALL MPI_send(my_n_orb_u_Bi,pksz_Bi,MPI_double_precision,0,  &
                 5000+myid,MPI_comm_world,ierr)
   CALL MPI_send(my_n_orb_d_Bi,pksz_Bi,MPI_double_precision,0,  &
                 6000+myid,MPI_comm_world,ierr)
   CALL MPI_send(my_n_orb_u,pksz,MPI_double_precision,0,  &
                 1000+myid,MPI_comm_world,ierr)
   CALL MPI_send(my_n_orb_d,pksz,MPI_double_precision,0,  &
                 2000+myid,MPI_comm_world,ierr)
   CALL MPI_send(my_Sx,1,MPI_double_precision,0,       &
                   3000+myid,MPI_comm_world,ierr)
   CALL MPI_send(my_Sy,1,MPI_double_precision,0,       &
                   4000+myid,MPI_comm_world,ierr)

  END IF


  CALL MPI_bcast(n_orb_u_Bi,pksz_Bi,MPI_double_precision,0,MPI_COMM_WORLD,ierr)
  CALL MPI_bcast(n_orb_d_Bi,pksz_Bi,MPI_double_precision,0,MPI_COMM_WORLD,ierr)
  CALL MPI_bcast(n_orb_u,pksz,MPI_double_precision,0,MPI_COMM_WORLD,ierr)
  CALL MPI_bcast(n_orb_d,pksz,MPI_double_precision,0,MPI_COMM_WORLD,ierr)
  CALL MPI_bcast(Sx,Nadatoms,MPI_double_precision,0,MPI_COMM_WORLD,ierr)
  CALL MPI_bcast(Sy,Nadatoms,MPI_double_precision,0,MPI_COMM_WORLD,ierr)

  n_orb_u_Bi = n_orb_u_Bi/pi + .5d0
  n_orb_d_Bi = n_orb_d_Bi/pi + .5d0
  n_orb_u = n_orb_u/pi + .5d0
  n_orb_d = n_orb_d/pi + .5d0
  Sx = Sx/pi
  Sy = Sy/pi



   n_t = SUM( n_orb_u + n_orb_d )
   DO is=1,Nadatoms
   m_f(is) = SUM(n_orb_u(is,Norb_ad-4:Norb_ad) - n_orb_d(is,Norb_ad-4:Norb_ad))
   END DO

  !fvec(1:Nadatoms) =  xnofel - n_t
  !fvec(1:Nadatoms) =  xnofel - SUM(n_orb_u(1,2:10)+n_orb_d(1,2:10))
  !fvec(Nadatoms+1:2*Nadatoms) = m_i - m_f
  !fvec(2*Nadatoms+1:3*Nadatoms) = Sx
  !fvec(3*Nadatoms+1:4*Nadatoms) = Sy

  n_up_Bi = 0.d0
  n_dn_Bi = 0.d0
  DO is=1,N_Bi
  n_up_Bi = n_up_Bi + SUM(n_orb_u_Bi(is,:)) 
  n_dn_Bi = n_dn_Bi + SUM(n_orb_d_Bi(is,:)) 
  END DO
  n_T_Bi = n_up_Bi+n_dn_Bi

  IF (myid==0) THEN
   WRITE(*,"('------------------')")
   WRITE(*,*)"Occupancies per orbital; eF = ",ef
   WRITE(*,*)
   WRITE(*,*)n_t,m_f
   WRITE(*,*)
   WRITE(*,*)"Sx,Sy"
   WRITE(*,*)Sx,Sy
   WRITE(*,*)
   WRITE(*,*)
   DO is=1,Nadatoms
   WRITE(*,*)
   WRITE(*,"(' s : ',f22.11)")n_orb_u(is,1)+n_orb_d(is,1)
   WRITE(*,"(' p : ',f22.11)")SUM(n_orb_u(is,2:4)+n_orb_d(is,2:4))
   WRITE(*,"(' d : ',f22.11)")SUM(n_orb_u(is,5:9)+n_orb_d(is,5:9))
   WRITE(*,"('d1 : ',3(f22.11,1x))")n_orb_u(is,5),n_orb_d(is,5),(n_orb_u(is,5)+n_orb_d(is,5))
   WRITE(*,"('d2 : ',3(f22.11,1x))")n_orb_u(is,6),n_orb_d(is,6),(n_orb_u(is,6)+n_orb_d(is,6))
   WRITE(*,"('d3 : ',3(f22.11,1x))")n_orb_u(is,7),n_orb_d(is,7),(n_orb_u(is,7)+n_orb_d(is,7))
   WRITE(*,"('d4 : ',3(f22.11,1x))")n_orb_u(is,8),n_orb_d(is,8),(n_orb_u(is,8)+n_orb_d(is,8))
   WRITE(*,"('d5 : ',3(f22.11,1x))")n_orb_u(is,9),n_orb_d(is,9),(n_orb_u(is,9)+n_orb_d(is,9))
   WRITE(*,"('d polarization: ',f22.11)")SUM(n_orb_u(is,5:9)-n_orb_d(is,5:9))
   END DO
   WRITE(*,*)
   WRITE(*,*)
   DO is=1,N_Bi
   WRITE(*,"('Bi ',I0,' occupancy: ',3(f18.6,1x))")is, &
   SUM(n_orb_u_Bi(is,:)+n_orb_d_Bi(is,:))
   WRITE(*,*)
   END DO
  END IF

  RETURN


END SUBROUTINE fcn





SUBROUTINE sumC_linear(ec,gout_Bi,gout_TM)

 USE f90_kind
 USE constants_u
 USE lattice
 USE sum_k


 COMPLEX(double) :: ec
 COMPLEX(double), DIMENSION(N_Bi,N_Bi,2*Norb_Bi,2*Norb_Bi), &
                                            INTENT(OUT) :: gout_Bi
 COMPLEX(double), DIMENSION(Nadatoms,Nadatoms,2*Norb_ad,2*Norb_ad), &
                                            INTENT(OUT) :: gout_TM

 COMPLEX(double), DIMENSION(N_Bi,N_Bi,2*Norb_Bi,2*Norb_Bi) :: dG_Bi
 COMPLEX(double), DIMENSION(Nadatoms,Nadatoms,2*Norb_ad,2*Norb_ad) :: dG_TM

 INTEGER :: m,n,nkmax,iz
 REAL(double) :: beta
 REAL(double) :: auxk,auxw,wk,aukx,auky,akx,aky

 REAL(double), DIMENSION(2) :: qv

 ! Generation of the Cunningham points for a primitive rectangular lattice
 ! They are generated in units of 2pi/a, where a is the fcc lattice constant.
 
  nkmax = ncp

  auxk  = DBLE( 2*ncp )
  auxw  = DBLE( nkmax )
      
  gout_Bi = zero
  gout_TM = zero

  !$OMP PARALLEL
  !$OMP DO PRIVATE(aukx, qv, dG_Bi, dG_TM) REDUCTION(+:gout_Bi,gout_TM)
  DO n = 1, nkmax
  aukx = 0.5d0*DBLE(2*n-1)/auxk

    qv = [aukx,0.d0]
          
    CALL green(ec,qv,dG_Bi,dG_TM)

    gout_Bi = gout_Bi + dG_Bi/auxw/2d0
    gout_TM = gout_TM + dG_TM/auxw/2d0

    qv = [-aukx,0.d0]
          
    CALL green(ec,qv,dG_Bi,dG_TM)

    gout_Bi = gout_Bi + dG_Bi/auxw/2d0
    gout_TM = gout_TM + dG_TM/auxw/2d0


  END DO
  !$OMP END DO
  !$OMP END PARALLEL

  RETURN

END SUBROUTINE sumC_linear



