SUBROUTINE invers(A,N)

   USE f90_kind
   USE mkl95_lapack

   INTEGER :: N
   INTEGER, DIMENSION(N) :: ipiv
   INTEGER :: info
   COMPLEX(double), DIMENSION(N,N) :: A

   CALL getrf(A,ipiv,info)
   IF (info/=0) WRITE(1200,*)"GETRF info",info
   info=0
   CALL getri(A,ipiv,info)
   IF (info/=0) WRITE(1201,*)"GETRI info",info

   

   RETURN

END SUBROUTINE invers
