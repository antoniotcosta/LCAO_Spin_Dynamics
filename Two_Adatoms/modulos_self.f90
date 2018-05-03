MODULE  lattice

 USE f90_kind

 ! Nadatoms is the number of adatom in a unit cell
 INTEGER, PARAMETER :: N_Bi = 108, Nadatoms = 1, N_uc = N_Bi+NAdatoms+6 
 INTEGER, PARAMETER :: Nadtotal = 2 ! number of adatoms in the system
 INTEGER, PARAMETER :: Norb_Bi = 9, Norb_ad = 9
 INTEGER, PARAMETER :: dimh = N_Bi*Norb_Bi + Nadatoms*Norb_ad + 6
 INTEGER, PARAMETER :: dimh_p1 = N_Bi*Norb_Bi/3 + 2
 INTEGER, PARAMETER :: dimh_p3 = N_Bi*Norb_Bi + 6
 INTEGER, PARAMETER :: dimG = 2*dimH, dimGp3 = 2*dimH_p3
 INTEGER, PARAMETER :: nnmax=8
 INTEGER :: dimHBi
 REAL(double), DIMENSION(3) :: b1,b2
 REAL(double), DIMENSION(N_uc,3) :: r0

 INTEGER :: numbneighb ! Read in routine hamiltonian

 REAL(double), DIMENSION(nnmax,3) :: dij

 INTEGER :: Nspacer

END MODULE lattice



MODULE the_hamiltonian

 USE lattice

 COMPLEX(double), DIMENSION(dimH_p3,dimH_p3) :: H00_pristine,H01_pristine
 COMPLEX(double), DIMENSION(nnmax,dimH_p1,dimH_p1) :: Hij
 COMPLEX(double), DIMENSION(dimH,dimH) :: H00_central
 COMPLEX(double), DIMENSION(2*dimH_p3,dimG) :: Hllc,Hrlc
 COMPLEX(double), DIMENSION(dimG,2*dimH_p3) :: Hcll,Hcrl
 COMPLEX(double), DIMENSION(2*dimH_p3,2*dimH_p3) :: mum_l,Hp01,Hp10
 COMPLEX(double), DIMENSION(dimG,dimG) :: mum,H00,HAB,HBA

 REAL(double) :: deltaEF_leads,deltaEF_central

END MODULE the_hamiltonian



module splitting

 USE lattice

 REAL(double), DIMENSION(Nadtotal) :: e0d
 COMPLEX(double), DIMENSION(Nadtotal,2,2) :: hdel_A
 REAL(double) :: hw0,U

end module splitting



MODULE ef_and_eta

 USE f90_kind

 REAL(double) :: ef,eta

END MODULE ef_and_eta



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



MODULE param_fcn

 USE lattice

 REAL(double), DIMENSION(Nadtotal) :: xnofel

END MODULE param_fcn
 

MODULE MPI_pars

 USE f90_kind

 INTEGER :: myid,numprocs
   
END MODULE MPI_pars



!--------------------------End of MODULES section--------------------------

