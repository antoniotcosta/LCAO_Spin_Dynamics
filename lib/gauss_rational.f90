SUBROUTINE gauss_rational(a,eta,x,w,nx)

 USE f90_kind
 USE constants_u

 REAL(double), INTENT(IN) :: a,eta
 INTEGER, INTENT(IN) :: nx
 REAL(double), DIMENSION(nx), INTENT(OUT) :: x,w

 INTEGER :: itype,ifail

 EXTERNAL d01bay

 itype = 0
 ifail = 0

 CALL d01bbf(d01bay,a,eta,itype,nx,w,x,ifail)

 IF (ifail/=0) THEN 
    WRITE(*,*)"problema na gauss_rational; ifail=",ifail
    STOP
 END IF

END SUBROUTINE gauss_rational
