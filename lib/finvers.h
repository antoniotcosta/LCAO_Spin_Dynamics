INTERFACE 
   FUNCTION finvers(matriz,dim)
      USE f90_kind

      INTEGER, INTENT(IN) :: dim
      COMPLEX(DOUBLE), DIMENSION(dim,dim), INTENT(IN) :: matriz
      COMPLEX(DOUBLE), DIMENSION(dim,dim) :: finvers
      
   END FUNCTION finvers

END INTERFACE 
