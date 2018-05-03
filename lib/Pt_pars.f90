SUBROUTINE Pt_pars

 USE f90_kind
 USE tight_binding
 USE Fermi_Energies


  s0(1)  = 0.77393d0 - EF_Pt
  p0(1)  = 1.48490d0 - EF_Pt
  d0e(1) = 0.45468d0 - EF_Pt
  d0t(1) = 0.43701d0 - EF_Pt

!     first n.n.

  sss(1,1) = -0.07835d0
  pps(1,1) =  0.18677d0
  ppp(1,1) = -0.02465d0
  dds(1,1) = -0.06856d0
  ddp(1,1) =  0.03528d0
  ddd(1,1) = -0.00588d0
  sps(1,1) =  0.11192d0
  sds(1,1) = -0.06197d0
  pds(1,1) = -0.08564d0
  pdp(1,1) =  0.02446d0

!     second n.n.

  sss(1,2) =  0.00157d0
  pps(1,2) =  0.03487d0
  ppp(1,2) = -0.00958d0
  dds(1,2) = -0.00300d0
  ddp(1,2) =  0.00226d0
  ddd(1,2) = -0.00029d0
  sps(1,2) =  0.00074d0
  sds(1,2) = -0.00770d0
  pds(1,2) = -0.01321d0
  pdp(1,2) =  0.00522d0

 RETURN

END SUBROUTINE Pt_pars
