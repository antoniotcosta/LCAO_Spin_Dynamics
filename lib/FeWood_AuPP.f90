SUBROUTINE FeWood_AuPP

  USE f90_kind
  USE tight_binding
  USE fermi

  !     two centre integrals; 
  !     Ag and Au (FCC): Papaconstantopoulos parameters. 
  !     Fe parameters fitting a three-centre band structure by Wood
  !
  !     The first index in the tight binding parameter arrays refers to
  !         1: Fe 
  !         2: Au

  ! Fe
  ! on site
  

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

  !   Au
  !   on site
  
  
  s0(2) =  0.56220d0 + Au_shift
  p0(2) =  1.27897d0 + Au_shift
  d0t(2) = 0.26097d0 + Au_shift
  d0e(2) = 0.25309d0 + Au_shift

  !  first n.n.

  sss(2,1) = -0.06680d0
  sps(2,1) =  0.09721d0
  pps(2,1) =  0.17866d0
  ppp(2,1) = -0.01645d0
  sds(2,1) = -0.04722d0
  pds(2,1) = -0.06399d0
  pdp(2,1) =  0.01896d0
  dds(2,1) = -0.04971d0
  ddp(2,1) =  0.02624d0
  ddd(2,1) = -0.00457d0

  !  second n.n.

  sss(2,2) =  0.00277d0
  sps(2,2) =  0.00261d0
  pps(2,2) =  0.03707d0
  ppp(2,2) = -0.01025d0
  sds(2,2) = -0.00784d0
  pds(2,2) = -0.00762d0
  pdp(2,2) =  0.00470d0
  dds(2,2) = -0.00305d0
  ddp(2,2) =  0.00240d0
  ddd(2,2) = -0.00057d0

  RETURN

END SUBROUTINE FeWood_AuPP
