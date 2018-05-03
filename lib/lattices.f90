MODULE lattice_fcc

  USE f90_kind

  INTEGER :: n0,n01,n02,n1,n2
  REAL(DOUBLE), DIMENSION(10,3) :: r0,c0,r1,c1,r2,c2

CONTAINS

  SUBROUTINE fcc001

    USE constants_u

    INTEGER :: i
    
    REAL(DOUBLE) :: length_r

    REAL(DOUBLE), DIMENSION(10,3) :: r0_r, r1_r

    !   in plane:

    n01 = 4
    n02 = 4

    r0(1,:) = [1.d0, 1.d0, 0.d0]
    r0(2,:) = [1.d0, -1.d0, 0.d0]
    r0(3,:) = [-1.d0, 1.d0, 0.d0]
    r0(4,:) = [-1.d0, -1.d0, 0.d0]

    r0(5,:) = [2.d0, 0.d0, 0.d0]
    r0(6,:) = [0.d0, 2.d0, 0.d0]
    r0(7,:) = [-2.d0, 0.d0, 0.d0]
    r0(8,:) = [0.d0, -2.d0, 0.d0]

    !  inter plane

    n1 = 4
    n2 = 1

    r1(1,:) = [1.d0, 0.d0, 1.d0]
    r1(2,:) = [0.d0, 1.d0, 1.d0]
    r1(3,:) = [0.d0, -1.d0, 1.d0]
    r1(4,:) = [-1.d0, 0.d0, 1.d0]

    r2(1,:) = [0.d0, 0.d0, 2.d0]

    DO i=1,n01+n02
       length_r = SQRT( DOT_PRODUCT(r0(i,:),r0(i,:)) )
       c0(i,:) = r0(i,:)/length_r
    END DO

    length_r = SQRT( DOT_PRODUCT(r1(1,:),r1(1,:)) )
    c1 = r1/length_r

    length_r = SQRT( DOT_PRODUCT(r2(1,:),r2(1,:)) )
    c2 = r2/length_r

    RETURN

  END SUBROUTINE fcc001


END MODULE lattice_fcc
