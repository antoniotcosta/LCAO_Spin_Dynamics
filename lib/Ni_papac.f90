subroutine Ni_papac

 use f90_kind
 use fermi_and_eta
 use tight_binding
 


! paramagnetic Ni (Papaconstantopoulos)

  s0(3)  = 1.15526d0 - ef_Ni + ef_Cu
  p0(3)  = 1.60587d0 - ef_Ni + ef_Cu
  d0(3) = 0.5d0*(0.56068d0 + 0.55753d0) - ef_Ni + ef_Cu
  !d0t(3) = 0.56068d0 - ef_Ni + ef_Cu
  !d0e(3) = 0.55753d0 - ef_Ni + ef_Cu

!     first n.n.
      
  sss(3,1) = -0.09525d0
  pps(3,1) =  0.21708d0
  ppp(3,1) =  0.01660d0
  dds(3,1) = -0.03712d0
  ddp(3,1) =  0.02629d0
  ddd(3,1) = -0.00600d0
  sps(3,1) =  0.14003d0
  sds(3,1) = -0.03880d0
  pds(3,1) = -0.04400d0
  pdp(3,1) =  0.02377d0
   
!     second n.n.

  sss(3,2) = -0.00065d0
  pps(3,2) =  0.06220d0
  ppp(3,2) =  0.00682d0
  dds(3,2) = -0.00651d0
  ddp(3,2) =  0.00344d0
  ddd(3,2) = -0.00027d0
  sps(3,2) =  0.01441d0
  sds(3,2) = -0.01015d0
  pds(3,2) = -0.01012d0
  pdp(3,2) =  0.00510d0

  return

end 
