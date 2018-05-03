MODULE  lattice

 USE f90_kind

 INTEGER, PARAMETER :: N_Bi = 16, Nadatoms = 1, N_uc = N_Bi+NAdatoms
 INTEGER, PARAMETER :: Nadtotal = 1 ! number of adatoms in the system
 INTEGER, PARAMETER :: Norb_Bi = 9, Norb_ad = 9
 INTEGER, PARAMETER :: dimh = N_Bi*Norb_Bi + Nadatoms*Norb_ad
 INTEGER, PARAMETER :: dimG = 2*dimH
 INTEGER, PARAMETER :: nnmax=8
 INTEGER :: dimHBi
 REAL(double), DIMENSION(3) :: b1,b2
 REAL(double), DIMENSION(N_uc,3) :: r0

 INTEGER :: numbneighb ! Read in routine hamiltonian

 REAL(double), DIMENSION(nnmax,3) :: dij

END MODULE lattice



MODULE the_hamiltonian

 USE lattice

 COMPLEX(double), DIMENSION(dimH,dimH) :: H00_central
 COMPLEX(double), DIMENSION(dimG,dimG) :: Hllc,Hrlc
 COMPLEX(double), DIMENSION(dimG,dimG) :: Hcll,Hcrl
 COMPLEX(double), DIMENSION(dimG,dimG) :: mum_l,Hp01,Hp10,tl,tld
 COMPLEX(double), DIMENSION(dimG,dimG) :: mum,H00,H01_central,H10_central

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



MODULE MPI_pars

 USE f90_kind

 INTEGER :: myid,numprocs
   
END MODULE MPI_pars

