SUBROUTINE bcc_bulk

 USE f90_kind
 USE constants_u
 USE lattice

 INTEGER :: i

 REAL(double) :: aux0

 n0 = 14

 r0(1,:) = [ 2.d0, 0.d0, 0.d0 ]

 r0(2,:) = [ 0.d0, 2.d0, 0.d0 ]

 r0(3,:) = [-2.d0, 0.d0, 0.d0 ]

 r0(4,:) = [ 0.d0,-2.d0, 0.d0 ]

 r0(5,:) = [ 0.d0, 0.d0, 2.d0 ]

 r0(6,:) = [ 0.d0, 0.d0,-2.d0 ]

 r0(7,:) = [ 1.d0, 1.d0, 1.d0 ] 

 r0(8,:) = [-1.d0,-1.d0,-1.d0 ]

 r0(9,:) = [-1.d0, 1.d0, 1.d0 ]

 r0(10,:) = [ 1.d0,-1.d0,-1.d0 ]

 r0(11,:) = [ 1.d0,-1.d0, 1.d0 ]
 
 r0(12,:) = [-1.d0,1.d0,-1.d0 ]
 
 r0(13,:) = [ 1.d0, 1.d0,-1.d0 ]
 
 r0(14,:) = [-1.d0,-1.d0, 1.d0 ]

 DO i=1,n0
    c0 = r0/SQRT(DOT_PRODUCT(r0(i,:),r0(i,:)))
 END DO

 RETURN

END SUBROUTINE bcc_bulk
