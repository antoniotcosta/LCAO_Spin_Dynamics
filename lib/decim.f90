subroutine decim( w, t, td, ws, wo, wb, cri, n, itmax)

 use f90_kind

 integer :: n,it,itmax
 real(double) :: cri,test
 complex(double), dimension(n,n) ::  w,t,td,wb,ws,wo,wi,twi,wit,witd,twitd,   &
                                     tdwit,tn,tdn,ts,tds

  ts = t    ! save h01 and h10
  tds = td  ! 

  wb = w
  ws = w
  wo = w

  do it=1,itmax                    

   wi = wb
   call invers(wi,n)
      
   twi  = matmul(t,wi)
   wit  = matmul(wi,t)            
   witd = matmul(wi,td)
   
   twitd  = matmul(twi,td)
   tdwit  = matmul(td,wit)

   ws   = ws - twitd 
   wo   = wo - tdwit      
   wb   = wb - twitd - tdwit
   tn   = matmul(t,wit)
   tdn  = matmul(td,witd) 
   
   t  = tn
   td = tdn
   test =  sum(abs(t))

   if(test < cri) then
    t = ts
    td = tds
    return 
   end if

  end do

  write(6,600) itmax, cri

  600 format(1x,' it was not possible to achieve convergency',/,  &
            '       itmax =',i3,'   cri = ',e10.3)    

  stop

end subroutine decim            
