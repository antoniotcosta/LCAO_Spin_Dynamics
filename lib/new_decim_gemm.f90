subroutine decim( w, ts, tds, ws, wo, wb, cri, n, itmax)

 use f90_kind

 INTEGER, INTENT(IN) :: n
 COMPLEX(double), DIMENSION(n,n), INTENT(IN) :: w,ts,tds
 COMPLEX(double), DIMENSION(n,n), INTENT(OUT) :: ws,wo,wb
 REAL(double), INTENT(IN) :: cri
 INTEGER, INTENT(IN) :: itmax

 integer :: it,mu
 real(double) :: test
 complex(double), dimension(n,n) ::  wi,twi,wit,witd,twitd,tdwit,tn,tdn,t,td
 
 CHARACTER(LEN=1) :: transa,transb,transc
 INTEGER :: na,nb,nc,lda,ldb,ldc
 COMPLEX(double) :: alpha,beta
 COMPLEX(double), DIMENSION(n,n) :: C


 na = n
 nb = n
 nc = n
 lda = n
 ldb = n
 ldc = n
 
 transa='n'
 transb='n'
 alpha = CMPLX(1d0,0d0,double)
 beta = CMPLX(0d0,0d0,double)
 C = CMPLX(0d0,0d0,double)

  t = ts    
  td = tds 

  wb = w
  ws = w
  wo = w

  do it=1,itmax                    

   wi = wb
   call invers(wi,n)
      
   ! twi  = matmul(t,wi)
   CALL zgemm(transa, transb, na, nb, nc, alpha, t, lda, wi, ldb, beta, twi, ldc)
   !wit  = matmul(wi,t)            
   CALL zgemm(transa, transb, na, nb, nc, alpha, wi, lda, t, ldb, beta, wit, ldc)
   !witd = matmul(wi,td)
   CALL zgemm(transa, transb, na, nb, nc, alpha, wi, lda, td, ldb, beta, witd, ldc)
   
   !twitd  = matmul(twi,td)
   CALL zgemm(transa, transb, na, nb, nc, alpha, twi, lda, td, ldb, beta, twitd, ldc) 
   !tdwit  = matmul(td,wit)
   CALL zgemm(transa, transb, na, nb, nc, alpha, td, lda, wit, ldb, beta, tdwit, ldc)

   ws   = ws - twitd 
   wo   = wo - tdwit      
   wb   = wb - twitd - tdwit
   !tn   = matmul(t,wit)
   CALL zgemm(transa, transb, na, nb, nc, alpha, t, lda, wit, ldb, beta, tn, ldc)
   !tdn  = matmul(td,witd) 
   CALL zgemm(transa, transb, na, nb, nc, alpha, td, lda, witd, ldb, beta, tdn, ldc)
   
   t  = tn
   td = tdn
   test =  sum(abs(t))

   if(test < cri) then
    return 
   end if

  end do

  write(*,"('it was not possible to achieve convergency',/,'itmax =',i3,'   cri = ',e10.3)")itmax,cri

  WRITE(*,"('Test = ',f21.11,1x,' eta = ',f21.11)")test,AIMAG(w(1,1))

  stop

end subroutine decim            
