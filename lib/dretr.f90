function dretr(gg)

 use f90_kind

 integer :: i
 real(double) :: suma,dretr
 complex(double) :: gg(9,9)

  suma = 0.0d0
  do i=1,9
   suma = suma + real(gg(i,i))
  end do 
  dretr = suma

  return

end function dretr   
