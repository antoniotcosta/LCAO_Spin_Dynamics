INCLUDE '/home/antc/lib/f90_kind.f90'
INCLUDE '/home/antc/lib/constants_u.f90'



MODULE  lattice

 USE f90_kind

 INTEGER, PARAMETER :: N_Bi = 108, Nadatoms = 1, N_uc = N_Bi+NAdatoms+6 
 INTEGER, PARAMETER :: Norb_Bi = 9, Norb_ad = 9
 INTEGER, PARAMETER :: dimh = N_Bi*Norb_Bi + Nadatoms*Norb_ad + 6
 INTEGER, PARAMETER :: dimh_p1 = N_Bi*Norb_Bi/3 + 2
 INTEGER, PARAMETER :: dimh_p3 = N_Bi*Norb_Bi + 6
 INTEGER, PARAMETER :: dimG = 2*dimH
 INTEGER, PARAMETER :: nnmax=8
 REAL(double), DIMENSION(3) :: b1,b2
 REAL(double), DIMENSION(N_uc,3) :: r0

 INTEGER :: numbneighb ! Read in routine hamiltonian

 REAL(double), DIMENSION(nnmax,3) :: dij

END MODULE lattice



MODULE the_hamiltonian

 USE lattice

 COMPLEX(double), DIMENSION(dimH_p3,dimH_p3) :: H00_pristine,H01_pristine
 COMPLEX(double), DIMENSION(nnmax,dimH_p1,dimH_p1) :: Hij
 COMPLEX(double), DIMENSION(dimH,dimH) :: H00_central
 COMPLEX(double), DIMENSION(2*dimH_p3,dimG) :: Hllc,Hrlc
 COMPLEX(double), DIMENSION(dimG,2*dimH_p3) :: Hcll,Hcrl
 COMPLEX(double), DIMENSION(2*dimH_p3,2*dimH_p3) :: mum_l,tl,tld
 COMPLEX(double), DIMENSION(dimG,dimG) :: mum

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



module splitting

 USE lattice

 REAL(double), DIMENSION(Nadatoms) :: hdel_A,e0d
 REAL(double) :: hw0,U

end module splitting


MODULE MPI_pars

 USE f90_kind

 INTEGER :: myid,numprocs
   
END MODULE MPI_pars



!--------------------------End of MODULES section--------------------------



PROGRAM DOS_graphene_with_N_adatoms

 USE lattice
 USE constants_u
 USE ef_and_eta
 USE sum_k
 USE lambda_SOC
 USE the_hamiltonian
 USE the_LS_matrix
 USE the_LS_p_matrix
 USE splitting
 USE MPI_pars
 USE MPI

 REAL(double) :: tol,xtol
 REAL(double) :: UeV,BZeeman

 INTEGER :: i,j,nw
 REAL(double) :: wmin,wmax,dw
 REAL(double), DIMENSION(Nadatoms) :: m_i_A

 CHARACTER(LEN=255) :: nomex,nomey,nomez

 REAL(double) :: theta0,phi0

 COMPLEX(double), DIMENSION(N_Bi,N_Bi,18,18) :: rhoij_Bi

 COMPLEX(double), DIMENSION(N_Bi,N_Bi,9,9) :: H00_Bi
 REAL(double), DIMENSION(18,19:108) :: Jsz,Jsx,Jsy
 INTEGER ::  ileft,iright
 INTEGER :: is,js,il,ir,mu,nu,mu0,muf,nu0,nuf,ief

 INTEGER :: IERR


  ! MPI Initialization
  CALL MPI_init(IERR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

  pi = 4.d0*ATAN(1.d0)
  sq2 = SQRT(2.d0)
  sq3 = SQRT(3.d0)

  OPEN(unit=99,file='spincurrent_entrada.in',status='old')

  READ(99,*)eta
  READ(99,*)UeV 
  U = UeV
  READ(99,*)BZeeman
  READ(99,*)e0d,m_i_A
  READ(99,*)theta0,phi0
  READ(99,*)lambda_Bi_p,lambda_Bi_d
  READ(99,*)lambda_ad
  READ(99,*)wmin,wmax,nw


  CLOSE(99)

  hdel_A = .5d0*U*m_i_A

  CALL hamiltonian_leads
  CALL hamiltonian_central
  CALL LS_matrix_p(theta0,phi0)

  WRITE(1110+myid,*)"finished setting up hamiltonians and LS_p"

  IF (myid==0) THEN
  DO is=1,18
  WRITE(nomex,"('Jsx_',I0,'.dat')")is
  WRITE(nomey,"('Jsy_',I0,'.dat')")is
  WRITE(nomez,"('Jsz_',I0,'.dat')")is
  OPEN(unit=124+is,file=nomex,status="unknown",position="append")
  OPEN(unit=224+is,file=nomey,status="unknown",position="append")
  OPEN(unit=324+is,file=nomez,status="unknown",position="append")
  END DO
  END IF

  dw = (wmax-wmin)/DBLE(nw)
  DO ief = 1,nw
  ef = wmin + DBLE(ief)*dw

  CALL correlationfunctions(rhoij_Bi)
  
  WRITE(1110+myid,*)"exited correlationfunctions"

  IF (myid==0) THEN

   DO i=1,N_Bi
   DO j=1,N_Bi

   mu0 = 9*Nadatoms + (i-1)*9 + 1; muf = mu0 + 8
   nu0 = 9*Nadatoms + (j-1)*9 + 1; nuf = nu0 + 8
   H00_Bi(i,j,:,:) = H00_central(mu0:muf,nu0:nuf)

   END DO
   END DO

   Jsx = 0d0
   Jsy = 0d0
   Jsz = 0d0


   DO is=1,18
   DO js=19,108
   DO mu=1,9
   DO nu=1,9
   Jsz(is,js) = Jsz(is,js) - 2.d0*H00_Bi(is,js,mu,nu)*	&
           ( AIMAG(rhoij_Bi(is,js,mu,nu)) - AIMAG(rhoij_Bi(is,js,mu+9,nu+9)) )
   Jsx(is,js) = Jsx(is,js) + H00_Bi(is,js,mu,nu)* &
           ( AIMAG(rhoij_Bi(is,js,mu,nu+9)) + AIMAG(rhoij_Bi(is,js,mu+9,nu)) ) 
   Jsy(is,js) = Jsy(is,js) + H00_Bi(is,js,mu,nu)* &
           ( REAL(rhoij_Bi(is,js,mu,nu+9)) - REAL(rhoij_Bi(is,js,mu+9,nu)) ) 
   END DO
   END DO
   END DO
   END DO

    DO is=1,18
    WRITE(124+is,"(92(f21.10,1x))")ef,Jsx(is,19:108),SUM(Jsx(is,19:108))
    WRITE(224+is,"(92(f21.10,1x))")ef,Jsy(is,19:108),SUM(Jsy(is,19:108))
    WRITE(324+is,"(92(f21.10,1x))")ef,Jsz(is,19:108),SUM(Jsz(is,19:108))
    END DO


  END IF
  END DO

  IF (myid==0) THEN
  DO is=1,18
  CLOSE(124+is)
  CLOSE(224+is)
  CLOSE(324+is)
  END DO
  END IF
  
 
  CALL MPI_Finalize(IERR)
  IF (IERR.ne.0) then
     WRITE(*,*)"stinky MPI finilize..."
  END IF

  STOP

END PROGRAM DOS_graphene_with_N_adatoms



SUBROUTINE correlationfunctions(rhoij_Bi)

 USE lattice
 USE constants_u
 USE ef_and_eta
 USE MPI
 USE MPI_pars

 COMPLEX(double), DIMENSION(N_Bi,N_Bi,18,18), INTENT(out) :: rhoij_Bi

 INTEGER :: n_o_nproc,resto,acum_n
 INTEGER, DIMENSION(0:99) :: npppr
 INTEGER :: ninte,i
 INTEGER, PARAMETER :: nptint = 64
 REAL(double) :: y,weoy2
 REAL(double) :: xe(nptint),we(nptint)
 REAL(double) :: my_xe(0:99,nptint),my_we(0:99,nptint)
 COMPLEX(double) :: ec
 COMPLEX(double), DIMENSION(N_Bi,N_Bi,18,18) :: &
                  gout_Bi,sumG_Bi,my_sumG_Bi,their_sumG_Bi


 INTEGER, DIMENSION(MPI_STATUS_SIZE) :: stat
 INTEGER :: ierr

 INTEGER :: pksz_Bi = N_Bi*N_Bi*18*18

 INTEGER :: il,is,js,mu,nu

 
  WRITE(1110+myid,*)"started correlfnctns"
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
 
  my_sumG_Bi = CMPLX(0d0,0d0,double)

  
 WRITE(1110+myid,*)"I have to calculate ",npppr(myid)," energy points"

  DO i=1,npppr(myid)
   y  = 1.d0 - my_xe(myid,i)
   ec = CMPLX(ef,(my_xe(myid,i)+eta)/y)
   weoy2 = (1.d0+eta)*my_we(myid,i)/(y*y)

   WRITE(1110+myid,*)"going into green"
   CALL green(ec,gout_Bi)
   WRITE(1110+myid,*)"out of green"

   my_sumG_Bi = my_sumG_Bi + weoy2*gout_Bi

    WRITE(1110+myid,*)"did my part",myid,i

  END DO

  ! Collect and add up
  IF (myid == 0) THEN
   sumG_Bi = my_sumG_Bi

   DO i=1,numprocs-1
    CALL MPI_Recv(their_sumG_Bi,pksz_Bi,MPI_double_complex, &
                  i,1000+i,MPI_comm_world,stat,ierr) 

    sumG_Bi = sumG_Bi + their_sumG_Bi
         
   END DO
  ELSE
   CALL MPI_send(my_sumG_Bi,pksz_Bi,MPI_double_complex,0,  &
                 1000+myid,MPI_comm_world,ierr)

    WRITE(1110+myid,*)"sent results"


  END IF

  IF (myid == 0) THEN

   DO is=1,N_Bi
   DO js=1,N_Bi
   DO mu=1,18
   DO nu=1,18
   IF((is==js).AND.(mu==nu))THEN 
   rhoij_Bi(is,js,mu,nu) = .5d0*( sumG_Bi(js,is,nu,mu) + &
                     CONJG(sumG_Bi(is,js,mu,nu)) )/pi + CMPLX(.5d0,0d0,double)
   ELSE 
   rhoij_Bi(is,js,mu,nu) = .5d0*( sumG_Bi(js,is,nu,mu) + &
                     CONJG(sumG_Bi(is,js,mu,nu)) )/pi
   END IF
   END DO
   END DO
   END DO
   END DO

  END IF


  RETURN


END SUBROUTINE correlationfunctions
