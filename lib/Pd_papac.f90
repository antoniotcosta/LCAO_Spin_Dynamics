subroutine Pd_papac

 use f90_kind
 use fermi
 use tight_binding
 

  s0(1)  = 0.94261d0
  p0(1)  = 1.36110d0
  d0e(1) = 0.37285d0
  d0t(1) = 0.36265d0

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

! paramagnetic Fe (Papaconstantopoulos)

  s0(2)  = 1.20177d0 + ef_Pd - ef_Fe
  p0(2)  = 1.87251d0 + ef_Pd - ef_Fe
  d0t(2) = 0.66437d0 + ef_Pd - ef_Fe
  d0e(2) = 0.68817d0 + ef_Pd - ef_Fe

!     first n.n.
      
  sss(2,1) = -0.13944d0
  sps(2,1) =  0.17780d0
  pps(2,1) =  0.26810d0
  ppp(2,1) =  0.02971d0
  sds(2,1) = -0.06781d0
  pds(2,1) = -0.09308d0
  pdp(2,1) =  0.02089d0
  dds(2,1) = -0.05086d0
  ddp(2,1) =  0.03096d0
  ddd(2,1) = -0.00303d0
   
!     second n.n.

  sss(2,2) = -0.03141d0
  sps(2,2) =  0.07354d0
  pps(2,2) =  0.18848d0
  ppp(2,2) =  0.03907d0
  sds(2,2) = -0.03884d0
  pds(2,2) = -0.06028d0
  pdp(2,2) = -0.00383d0
  dds(2,2) = -0.03125d0
  ddp(2,2) =  0.00618d0
  ddd(2,2) =  0.00071d0

  return
end 
