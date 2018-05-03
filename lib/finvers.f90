FUNCTION finvers(matriz,dim)
 
 USE f90_kind

 INTEGER, INTENT(IN) :: dim
 COMPLEX(double), DIMENSION(dim,dim), INTENT(IN) :: matriz
 
 COMPLEX(double), DIMENSION(dim,dim) :: matriz_out,finvers

 INTEGER, DIMENSION(dim) :: IPIV
 INTEGER :: INFO
 INTEGER :: LWORK 

 COMPLEX(double), DIMENSION(4*dim) :: WORK


  matriz_out = matriz

  INFO = 0
  LWORK = dim*4
  CALL F07ARF(dim,dim,matriz_out,dim,IPIV,INFO)          
  IF (INFO/=0) then
   WRITE(*,*)'IFAIL=',INFO
  END IF
  CALL F07AWF(dim,matriz_out,dim,IPIV,WORK,LWORK,INFO)

  IF (INFO/=0) THEN
    WRITE(*,*)'IFAIL=',INFO
    STOP 'FUDEU!!! INVERSAO DA NAG FOI PRO CARALHO!'
  END IF

  finvers = matriz_out

  RETURN

END FUNCTION finvers
