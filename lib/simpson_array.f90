program trapezium_array

 use f90_kind

 integer :: i,it,nt,nump
 real(double) :: x0(10000)
 complex(double) :: y0(10000),result,kern0,kern1,Iz

 character(len=255) :: nome_in,nome_out

 nome_in = 'chi_12.dat"
 nome_out = 'chi_12_t.dat"

 Iz = cmplx(0.d0,1.d0,double)

 open(unit=11,file=nome_in,status='old')
 open(unit=21,file=nome_in,status='old')

 do i=1,10000

  read(11,*,end=666)x0(i),y0(i)

 end do

 close(11)
 
 666 nump = i
 
 t0 = 0.d0
 tf = 10.d0
 nt = 1000
 tstep = (tf-t0)/float(nt)
 
 do it = 0,nt
  t = t0+tstep*float(it)
  result = cmplx(0.d0,0.d0,double)
  do i=1,nump-1
  
   kern0 = exp(Iz*x0(i)*t)
   kern1 = exp(Iz*x0(i+1)*t)

   result = result + (x0(i+1)-x0(i))*(y0(i+1)*kern1+y0(i)*kern0)*0.5d0

  end do

  write(21,*)t,abs(result)**2.D0

 end do

 close(21)

 stop

end program trapezium_array
 
