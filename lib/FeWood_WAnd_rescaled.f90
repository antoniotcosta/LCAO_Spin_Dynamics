SUBROUTINE Fe_W_andpar

 USE f90_kind
 USE fermi
 USE tight_binding
 USE distance_scaling

 REAL(double) :: shft,cs,cp,cd,ds,dp,dd,ds2,dsp,dsd,dp2,dpd,dd2

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

  ! on site
      
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

  ! first n.n.
      
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
   
  ! second n.n.

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
  !  Aligning the Fermi energy of paramagnetic Fe with the value of Ef_W

  shft = ef_W - ef_Fe

  ! on site
      
  s0(2)  = 1.1734246637016703d0 + shft
  p0(2)  = 1.5466607694715542d0 + shft
  d0t(2) = 0.7835466941723471d0 + shft     
  d0e(2) = 0.7600964802524716d0 + shft

  ! rescaling for change in a_0 
  !  first n.n.
      
  sss(2,1) = -0.1007863439107165d0/( (1.d0 + delta_a) )
  sps(2,1) = -0.1483480847871167d0/( (1.d0 + delta_a)**2.d0 )
  pps(2,1) =  0.1955800299799059d0/( (1.d0 + delta_a)**3.d0 )
  ppp(2,1) = -2.5318074485867866d-02/( (1.d0 + delta_a)**3.d0 )
  sds(2,1) = -6.7483442463254467d-02/( (1.d0 + delta_a)**3.d0 )
  pds(2,1) =  7.9400304909909458d-02/( (1.d0 + delta_a)**4.d0 )
  pdp(2,1) = -2.6623049988494171d-02/( (1.d0 + delta_a)**4.d0 )
  dds(2,1) = -4.8755321112568868d-02/( (1.d0 + delta_a)**5.d0 )
  ddp(2,1) =  3.3120152475773490d-02/( (1.d0 + delta_a)**5.d0 )
  ddd(2,1) = -4.4762536733048697d-03/( (1.d0 + delta_a)**5.d0 )
  
  ! second n.n.

  sss(2,2) = -3.6180669066274485d-02/( (1.d0 + delta_a) )
  sps(2,2) = -5.3667311183090011d-02/( (1.d0 + delta_a)**2.d0 )
  pps(2,2) =  9.7099646969018613d-02/( (1.d0 + delta_a)**3.d0 )
  ppp(2,2) =  1.7373459086569659d-02/( (1.d0 + delta_a)**3.d0 )
  sds(2,2) = -3.4902759086422826d-02/( (1.d0 + delta_a)**3.d0 )
  pds(2,2) =  4.9187216447024916d-02/( (1.d0 + delta_a)**4.d0 )
  pdp(2,2) = -2.1618689500718141d-03/( (1.d0 + delta_a)**4.d0 )
  dds(2,2) = -3.0802113741360822d-02/( (1.d0 + delta_a)**5.d0 )
  ddp(2,2) =  8.6494842901658053d-03/( (1.d0 + delta_a)**5.d0 )
  ddd(2,2) = -6.0886323666680103d-04/( (1.d0 + delta_a)**5.d0 )
 
  ! interplane hopping should not be changed
  ! preserved the unscaled parameters in a 3rd position

  sss(3,1) = -0.1007863439107165d0
  sps(3,1) = -0.1483480847871167d0
  pps(3,1) =  0.1955800299799059d0
  ppp(3,1) = -2.5318074485867866d-02
  sds(3,1) = -6.7483442463254467d-02
  pds(3,1) =  7.9400304909909458d-02
  pdp(3,1) = -2.6623049988494171d-02
  dds(3,1) = -4.8755321112568868d-02
  ddp(3,1) =  3.3120152475773490d-02
  ddd(3,1) = -4.4762536733048697d-03
  
  ! second n.n.

  sss(3,2) = -3.6180669066274485d-02
  sps(3,2) = -5.3667311183090011d-02
  pps(3,2) =  9.7099646969018613d-02
  ppp(3,2) =  1.7373459086569659d-02
  sds(3,2) = -3.4902759086422826d-02
  pds(3,2) =  4.9187216447024916d-02
  pdp(3,2) = -2.1618689500718141d-03
  dds(3,2) = -3.0802113741360822d-02
  ddp(3,2) =  8.6494842901658053d-03
  ddd(3,2) = -6.0886323666680103d-04

  
  RETURN

END SUBROUTINE Fe_W_andpar
