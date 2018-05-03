subroutine gauleg(x1,x2,x,w,N)

 use f90_kind

 integer :: N,IFAIL,ITYPE
 real(double) :: x1,x2,C,D,x(N),w(N),ABCISS(N),WEIGHT(N)

 IFAIL = 0
 ITYPE = 0
 C=0.D0 ! not used
 D=0.D0 ! not used
 call D01BCF(ITYPE,x1,x2,C,D,N,WEIGHT,ABCISS,IFAIL)

 x = ABCISS
 w = WEIGHT

 return

end subroutine gauleg
