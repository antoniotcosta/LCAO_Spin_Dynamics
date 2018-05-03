SUBROUTINE gauss_laguerre(A,x,w,N)

 use f90_kind

 integer, intent(in) :: N
 real(double), intent(in) :: A
 real(double), dimension(N), intent(out) :: x,w

 integer :: IFAIL,ITYPE
 real(double) :: B,C,D
 real(double), dimension(N) :: ABCISS,WEIGHT

 IFAIL = 0
 ITYPE = -3  ! Gauss-Laguerre, adjusted weights
 B = 10.d0    ! exp(-B*x)
 C = 0.D0    ! |x-a|^C
 D=0.D0 ! not used

 call D01BCF(ITYPE,A,B,C,D,N,WEIGHT,ABCISS,IFAIL)

 x = ABCISS
 w = WEIGHT

 RETURN

end SUBROUTINE gauss_laguerre
