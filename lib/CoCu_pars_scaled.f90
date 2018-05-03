SUBROUTINE CoCu_pars

 USE f90_kind
 USE tight_binding
 USE distance_scaling

 REAL(double) :: cshift

!     two centre integrals; slater-koster parameters for spin polarized
!     fcc Cobalt up, Cobalt down and Copper spin bands. 
!     For Cobalt, the exchange field was determined self-consistently 
!     starting from paramagnetic fcc cobalt parameters given by 
!     Papaconstantopoulos.
!     We add a constant shift to s0, p0, and d0 so that the Fermi level
!     of Co coincides with that of Cu.
!
!     The first index in the tight binding parameter arrays refers to
!         1: Co up
!         2: Co down
!         3: Cu
!
!     Co up:
  
  cshift = .575530d0 - .715751d0
  s0(1) =  1.12946d0 + cshift
  p0(1) =  1.75262d0 + cshift
  d0(1) =  0.5d0*(0.60547d0 + 0.60445d0) + cshift 

!     first n.n.
  
  sss(1,1) = -0.09043d0
  sps(1,1) =  0.13649d0
  pps(1,1) =  0.23748d0
  ppp(1,1) = -0.00142d0
  sds(1,1) = -0.03806d0
  pds(1,1) = -0.04069d0
  pdp(1,1) =  0.02797d0
  dds(1,1) = -0.04213d0
  ddp(1,1) =  0.02976d0
  ddd(1,1) = -0.00684d0
      
!     second n.n.

  sss(1,2) = -0.00337d0
  sps(1,2) =  0.00135d0
  pps(1,2) =  0.02849d0
  ppp(1,2) =  0.01099d0
  sds(1,2) = -0.01119d0
  pds(1,2) = -0.01061d0
  pdp(1,2) =  0.01134d0
  dds(1,2) = -0.00759d0
  ddp(1,2) =  0.00495d0
  ddd(1,2) = -0.00016d0
      
! Scaled for distance variation
!     first n.n.
  
  sss(3,1) = -0.09043d0/( (1.d0 + delta_a) )
  sps(3,1) =  0.13649d0/( (1.d0 + delta_a)**4.d0 )
  pps(3,1) =  0.23748d0/( (1.d0 + delta_a)**3.d0 )
  ppp(3,1) = -0.00142d0/( (1.d0 + delta_a)**3.d0 )
  sds(3,1) = -0.03806d0/( (1.d0 + delta_a)**3.d0 )
  pds(3,1) = -0.04069d0/( (1.d0 + delta_a)**4.d0 )
  pdp(3,1) =  0.02797d0/( (1.d0 + delta_a)**4.d0 )
  dds(3,1) = -0.04213d0/( (1.d0 + delta_a)**5.d0 )
  ddp(3,1) =  0.02976d0/( (1.d0 + delta_a)**5.d0 )
  ddd(3,1) = -0.00684d0/( (1.d0 + delta_a)**5.d0 )
      
!     second n.n.

  sss(3,2) = -0.00337d0/( (1.d0 + delta_a) )
  sps(3,2) =  0.00135d0/( (1.d0 + delta_a)**2.d0 )
  pps(3,2) =  0.02849d0/( (1.d0 + delta_a)**3.d0 )
  ppp(3,2) =  0.01099d0/( (1.d0 + delta_a)**3.d0 )
  sds(3,2) = -0.01119d0/( (1.d0 + delta_a)**3.d0 )
  pds(3,2) = -0.01061d0/( (1.d0 + delta_a)**4.d0 )
  pdp(3,2) =  0.01134d0/( (1.d0 + delta_a)**4.d0 )
  dds(3,2) = -0.00759d0/( (1.d0 + delta_a)**5.d0 )
  ddp(3,2) =  0.00495d0/( (1.d0 + delta_a)**5.d0 )
  ddd(3,2) = -0.00016d0/( (1.d0 + delta_a)**5.d0 )

!     Cu:
  
  s0(2) =  0.79466d0
  p0(2) =  1.35351d0
  d0(2) =  0.5d0*(0.37307d0 + 0.37180d0)

!     first n.n.
      
  sss(2,1) = -0.07518d0
  sps(2,1) =  0.11571d0
  pps(2,1) =  0.19669d0
  ppp(2,1) =  0.01940d0
  sds(2,1) = -0.03107d0
  pds(2,1) = -0.03289d0
  pdp(2,1) =  0.01753d0
  dds(2,1) = -0.02566d0
  ddp(2,1) =  0.01800d0
  ddd(2,1) = -0.00408d0
      
!     second n.n.

  sss(2,2) = -0.00092d0
  sps(2,2) =  0.01221d0
  pps(2,2) =  0.05389d0
  ppp(2,2) =  0.00846d0
  sds(2,2) = -0.00852d0
  pds(2,2) = -0.00536d0
  pdp(2,2) =  0.00321d0
  dds(2,2) = -0.00451d0
  ddp(2,2) =  0.00241d0
  ddd(2,2) = -0.00029d0
      
  RETURN
END SUBROUTINE CoCu_pars
