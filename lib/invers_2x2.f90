SUBROUTINE invers2(A)

 USE f90_kind

 COMPLEX(double), DIMENSION(2,2) :: A,Atemp

 COMPLEX(double) :: det

 Atemp = A
 det = A(1,1)*A(2,2) - A(2,1)*A(1,2)

 A(1,1) =  Atemp(2,2)/det; A(1,2) = -Atemp(1,2)/det
 A(2,1) = -Atemp(2,1)/det; A(2,2) =  Atemp(1,1)/det

 RETURN

END SUBROUTINE invers2
