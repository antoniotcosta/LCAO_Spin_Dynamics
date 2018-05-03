subroutine decim( w, ts, tds, ws, wo, wb, cri, n, itmax)

 use f90_kind
 USE mpi_pars

 INTEGER, INTENT(IN) :: n
 COMPLEX(double), DIMENSION(n,n), INTENT(IN) :: w,ts,tds
 COMPLEX(double), DIMENSION(n,n), INTENT(OUT) :: ws,wo,wb
 REAL(double), INTENT(IN) :: cri
 INTEGER, INTENT(IN) :: itmax

 integer :: it,mu
 real(double) :: test
 complex(double), dimension(n,n) ::  wi,twi,wit,witd,twitd,tdwit,tn,tdn,t,td

  t = ts    
  td = tds 

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
    return 
   end if

  end do

  write(*,"('it was not possible to achieve convergency',/,'itmax =',i3,'   cri = ',e10.3)")itmax,cri

  WRITE(100+myid,*)"phodeu"
  DO mu=1,n
     WRITE(100+myid,*)AIMAG(w(mu,mu))
  END DO

  stop

end subroutine decim            
