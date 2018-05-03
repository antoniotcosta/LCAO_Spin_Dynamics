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

! TB parameters fitting a three-centre band structure by Wood

!     Aligning the Fermi energy of paramagnetic Fe with the value of Ef_W

  shft = ef_W - ef_Fe

      
   s0(2) = 1.2017709017  + shft
   p0(2) = 1.8725119829  + shft
   d0t(2) = 0.6881678104 + shft
   d0e(2) = 0.6643740535 + shft

   sss(2,1) =-0.1394413859
   pps(2,1) = 0.2681021988
   ppp(2,1) = 0.0297146384
   dds(2,1) =-0.0508569255
   ddp(2,1) = 0.0309574008
   ddd(2,1) =-0.0030320531
   sps(2,1) = 0.1777951121
   sds(2,1) =-0.0678095073
   pds(2,1) =-0.0930757448
   pdp(2,1) = 0.0208929181

   sss(2,2) =-0.0314096436
   pps(2,2) = 0.1884829849
   ppp(2,2) = 0.0390681326
   dds(2,2) =-0.0312470067
   ddp(2,2) = 0.0061819027
   ddd(2,2) = 0.0007075703
   sps(2,2) = 0.0735426247
   sds(2,2) =-0.0388437621
   pds(2,2) =-0.0602805056
   pdp(2,2) =-0.0038276755
 
  return

end subroutine Fe_W_andpar
