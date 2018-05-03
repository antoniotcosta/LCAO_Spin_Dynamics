subroutine Ni_Pd

 use f90_kind
 use fermi
 use tight_binding
 

  s0(1)  = 0.94261d0 - ef_Pd + ef_Ni
  p0(1)  = 1.36110d0 - ef_Pd + ef_Ni
  d0e(1) = 0.37285d0 - ef_Pd + ef_Ni
  d0t(1) = 0.36265d0 - ef_Pd + ef_Ni

!     first n.n.
      
  sss(1,1) = -0.07962d0
  sps(1,1) =  0.11332d0
  pps(1,1) =  0.17119d0  
  ppp(1,1) = -0.00540d0
  sds(1,1) = -0.04885d0
  pds(1,1) = -0.06563d0
  pdp(1,1) =  0.02124d0 
  dds(1,1) = -0.05216d0
  ddp(1,1) =  0.02878d0 
  ddd(1,1) = -0.00533d0
   
!     second n.n.

  sss(1,2) = -0.00105d0
  sps(1,2) =  0.01048d0
  pps(1,2) =  0.04282d0
  ppp(1,2) = -0.00044d0
  sds(1,2) = -0.00837d0
  pds(1,2) = -0.00738d0
  pdp(1,2) =  0.00351d0
  dds(1,2) = -0.00385d0
  ddp(1,2) =  0.00212d0
  ddd(1,2) = -0.00026d0

! paramagnetic Ni (Papaconstantopoulos)

  s0(2)  = 1.15526d0
  p0(2)  = 1.60587d0
  d0t(2) = 0.56068d0
  d0e(2) = 0.55753d0

!     first n.n.
      
  sss(2,1) = -0.09525d0
  pps(2,1) =  0.21708d0
  ppp(2,1) =  0.01660d0
  dds(2,1) = -0.03712d0
  ddp(2,1) =  0.02629d0
  ddd(2,1) = -0.00600d0
  sps(2,1) =  0.14003d0
  sds(2,1) = -0.03880d0
  pds(2,1) = -0.04400d0
  pdp(2,1) =  0.02377d0
   
!     second n.n.

  sss(2,2) = -0.00065d0
  pps(2,2) =  0.06220d0
  ppp(2,2) =  0.00682d0
  dds(2,2) = -0.00651d0
  ddp(2,2) =  0.00344d0
  ddd(2,2) = -0.00027d0
  sps(2,2) =  0.01441d0
  sds(2,2) = -0.01015d0
  pds(2,2) = -0.01012d0
  pdp(2,2) =  0.00510d0

  return

end 
