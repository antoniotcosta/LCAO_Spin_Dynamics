subroutine Fe_W_andpar

 use f90_kind
 use fermi
 use tight_binding

 real(double) :: shft,cs,cp,cd,ds,dp,dd,ds2,dsp,dsd,dp2,dpd,dd2

! TB parameters fitting a three-centre band structure by Wood

!     Aligning the Fermi energy of paramagnetic Fe with the value of Ef_W

!     on site
      
  s0(1)  = 1.1734246637016703d0
  p0(1)  = 1.5466607694715542d0
  d0t(1) = 0.7835466941723471d0
  d0e(1) = 0.7600964802524716d0

!     first n.n.
      
  sss(1,1) = -0.1007863439107165d0
  sps(1,1) = -0.1483480847871167d0
  pps(1,1) =  0.1955800299799059d0
  ppp(1,1) = -2.5318074485867866d-02
  sds(1,1) = -6.7483442463254467d-02
  pds(1,1) =  7.9400304909909458d-02
  pdp(1,1) = -2.6623049988494171d-02
  dds(1,1) = -4.8755321112568868d-02
  ddp(1,1) =  3.3120152475773490d-02
  ddd(1,1) = -4.4762536733048697d-03
  
!     second n.n.

  sss(1,2) = -3.6180669066274485d-02
  sps(1,2) = -5.3667311183090011d-02
  pps(1,2) =  9.7099646969018613d-02
  ppp(1,2) =  1.7373459086569659d-02
  sds(1,2) = -3.4902759086422826d-02
  pds(1,2) =  4.9187216447024916d-02
  pdp(1,2) = -2.1618689500718141d-03
  dds(1,2) = -3.0802113741360822d-02
  ddp(1,2) =  8.6494842901658053d-03
  ddd(1,2) = -6.0886323666680103d-04
 
  return

end subroutine Fe_W_andpar
