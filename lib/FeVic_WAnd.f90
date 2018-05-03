subroutine Fe_W_andpar

 use f90_kind
 use fermi
 use tight_binding

 real(double) :: shft,cs,cp,cd,ds,dp,dd,ds2,dsp,dsd,dp2,dpd,dd2

!     two centre integrals; Slater-Koster-Andersen parameters for  
!     BCC structure       
!
!     W:
!

  cs = -.329020d0
  cp =  .397044d0
  cd =  .023321d0
  ds =  .387912d0
  dp =  .232742d0
  dd =  .180593d0

!     on site
      
  ds2 = ds*ds
  dsp = ds*dp
  dsd = ds*dd
  dp2 = dp*dp
  dpd = dp*dd
  dd2 = dd*dd 

  s0(1)  = cs + ds2*3.09d0
  p0(1)  = cp + dp2*2.79d0
  d0t(1) = cd + dd2*2.71d0
  d0e(1) = cd + dd2*1.30d0

!     first n.n.
      
  sss(1,1) = -ds2*0.593d0
  sps(1,1) =  dsp*1.18d0
  pps(1,1) =  dp2*2.36d0
  ppp(1,1) = -dp2*0.36d0
  sds(1,1) = -dsd*1.42d0
  pds(1,1) = -dpd*2.93d0
  pdp(1,1) =  dpd*0.82d0
  dds(1,1) = -dd2*3.84d0
  ddp(1,1) =  dd2*1.85d0
  ddd(1,1) = -dd2*0.19d0
   
!     second n.n.

  sss(1,2) = -ds2*0.203d0
  sps(1,2) =  dsp*0.44d0
  pps(1,2) =  dp2*0.93d0
  ppp(1,2) = -dp2*0.05d0
  sds(1,2) = -dsd*0.60d0
  pds(1,2) = -dpd*1.29d0
  pdp(1,2) =  dpd*0.13d0
  dds(1,2) = -dd2*1.76d0
  ddp(1,2) =  dd2*0.36d0
  ddd(1,2) = -dd2*0.02d0

!     two centre integrals; slater-koster parameters for bcc
!     paramagnetic Iron from R H Victora "Magnetic and Electronic
!     Properties of Transition Metals and Overlayers"

 
!     Aligning the Fermi energy of paramagnetic Fe with the value of
!     Ef_W shft = ef_W - ef_Fe should be added to
!     s0(2),p0(2),d0t(2),d0e(2) on site

      shft = ef_W - ef_Fe
      
      s0(2)  =  0.514d0 + shft
      p0(2)  =  1.118d0 + shft
      d0t(2) = -0.089d0 + shft
      d0e(2) = -0.068d0 + shft

!     first n.n.
      
      sss(2,1) = -0.10605d0
      sps(2,1) = -0.19326d0
      pps(2,1) =  0.26932d0
      ppp(2,1) = -0.00672d0
      sds(2,1) = -0.08018d0
      pds(2,1) =  0.09915d0
      pdp(2,1) = -0.02467d0
      dds(2,1) = -0.05669d0
      ddp(2,1) =  0.03738d0
      ddd(2,1) = -0.00669d0
      
!     second n.n.

      sss(2,2) = -0.05276d0
      sps(2,2) = -0.06725d0
      pps(2,2) =  0.17199d0
      ppp(2,2) =  0.03242d0
      sds(2,2) = -0.03035d0
      pds(2,2) =  0.05202d0
      pdp(2,2) = -0.00092d0
      dds(2,2) = -0.03367d0
      ddp(2,2) =  0.01000d0
      ddd(2,2) = -0.00075d0

 
  return

end subroutine Fe_W_andpar
