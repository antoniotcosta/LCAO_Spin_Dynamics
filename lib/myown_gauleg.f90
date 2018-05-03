SUBROUTINE gauleg(x1,x2,x,w,n)
 USE f90_kind
 USE constants_u
     
      INTEGER, INTENT(IN) :: n 
      REAL(double), INTENT(IN) :: x1,x2
      REAL(double), DIMENSION(n), INTENT(OUT) :: x,w

      REAL(double), PARAMETER :: EPS = 3.d-14
      INTEGER :: i,j,m
      REAL(double) :: p1,p2,p3,pp,xl,xm,z,z1

      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      DO i=1,m
        z=cos(pi*(i-.25d0)/(n+.5d0))
        DO WHILE (ABS(z-z1).gt.EPS)
           p1=1.d0
           p2=0.d0
           DO j=1,n
              p3=p2
              p2=p1
              p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
           END DO
           pp=n*(z*p1-p2)/(z*z-1.d0)
           z1=z
           z=z1-p1/pp
        END DO 
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
      END DO

      RETURN

END SUBROUTINE gauleg

