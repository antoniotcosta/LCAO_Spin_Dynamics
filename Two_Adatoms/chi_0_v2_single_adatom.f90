SUBROUTINE dchi_HF_1(w,csi,chi0m)

 USE f90_kind
 USE constants_u
 USE ef_and_eta
 USE lattice

  REAL(double), INTENT(IN) :: w,csi
  COMPLEX(double), DIMENSION(4*5*5*Nadtotal,5*5*Nadtotal), INTENT(OUT) :: chi0m

  COMPLEX(double), DIMENSION(3,Nadtotal,Nadtotal,18,18) :: gout
  COMPLEX(double), DIMENSION(Nadtotal,Nadtotal,18,18) :: gtemp
  COMPLEX(double), DIMENSION(3,Nadtotal,Nadtotal,9,9) :: g_upup,g_updn,g_dnup,g_dndn

  COMPLEX(double) :: wc

  INTEGER :: a,b,ll,l,m,n,mu,mul,nu,nul

  ! csi: integration variable

  wc = CMPLX(ef+w,csi)
  CALL green(wc,gtemp)
  gout(1,:,:,:,:) = gtemp

  wc = CMPLX(ef,csi)
  CALL green(wc,gtemp)
  gout(2,:,:,:,:) = gtemp
  
  wc = CMPLX(ef-w,csi)
  CALL green(wc,gtemp)
  gout(3,:,:,:,:) = gtemp

  DO l=1,3
  DO is=1,Nadtotal
  DO js=1,Nadtotal
  g_upup(l,is,js,:,:) = gout(l,is,js,1:9,1:9)
  g_updn(l,is,js,:,:) = gout(l,is,js,1:9,10:18)
  g_dnup(l,is,js,:,:) = gout(l,is,js,10:18,1:9)
  g_dndn(l,is,js,:,:) = gout(l,is,js,10:18,10:18)
  END DO
  END DO
  END DO

  l = 1
  DO is=1,Nadtotal
  DO js=1,Nadtotal
  DO mu = 1,5
  DO mul = 1,5
  m = (l-1)*25*Nadtotal + (is-1)*25 + (mu-1)*5 + mul 
  DO nu = 1,5
  DO nul = 1,5
  n = (js-1)*25 + (nu-1)*5 + nul

  chi0m(m,n) = g_dndn(1,is,js,nu+4,mul+4)*g_upup(2,js,is,nul+4,mu+4) + &
               CONJG(g_upup(3,is,js,mu+4,nul+4)*g_dndn(2,js,is,mul+4,nu+4))
  END DO
  END DO
  END DO
  END DO
  END DO
  END DO

  l = 2
  DO is=1,Nadtotal
  DO js=1,Nadtotal
  DO mu = 1,5
  DO mul = 1,5
  m = (l-1)*25*Nadtotal + (is-1)*25 + (mu-1)*5 + mul 
  DO nu = 1,5
  DO nul = 1,5
  n = (js-1)*25 + (nu-1)*5 + nul

  chi0m(m,n) = g_updn(1,is,js,nu+4,mul+4)*g_upup(2,js,is,nul+4,mu+4) +  &
               CONJG(g_upup(3,is,js,mu+4,nul+4)*g_dnup(2,js,is,mul+4,nu+4))
  END DO
  END DO
  END DO
  END DO
  END DO
  END DO

  l = 3
  DO is=1,Nadtotal
  DO js=1,Nadtotal
  DO mu = 1,5
  DO mul = 1,5
  m = (l-1)*25*Nadtotal + (is-1)*25 + (mu-1)*5 + mul 
  DO nu = 1,5
  DO nul = 1,5
  n = (js-1)*25 + (nu-1)*5 + nul

  chi0m(m,n) = g_dndn(1,is,js,nu+4,mul+4)*g_updn(2,js,is,nul+4,mu+4) + &
              CONJG(g_dnup(3,is,js,mu+4,nul+4)*g_dndn(2,js,is,mul+4,nu+4))
  END DO
  END DO
  END DO
  END DO
  END DO
  END DO
 
  l = 4
  DO is=1,Nadtotal
  DO js=1,Nadtotal
  DO mu = 1,5
  DO mul = 1,5
  m = (l-1)*25*Nadtotal + (is-1)*25 + (mu-1)*5 + mul 
  DO nu = 1,5
  DO nul = 1,5
  n = (js-1)*25 + (nu-1)*5 + nul

  chi0m(m,n) = g_updn(1,is,js,nu+4,mul+4)*g_updn(2,js,is,nul+4,mu+4) + &
              CONJG(g_dnup(3,is,js,mu+4,nul+4)*g_dnup(2,js,is,mul+4,nu+4))
  END DO
  END DO
  END DO
  END DO
  END DO
  END DO

  chi0m = chi0m/(2.d0*pi)

  RETURN

END SUBROUTINE dchi_HF_1





SUBROUTINE dchi_HF_2(w,csi,chi0m)

 USE f90_kind
 USE constants_u
 USE ef_and_eta
 USE lattice


  REAL(double), INTENT(IN) :: w,csi
  COMPLEX(double), DIMENSION(4*5*5*Nadtotal,5*5*Nadtotal), INTENT(OUT) :: chi0m

  COMPLEX(double), DIMENSION(2,Nadtotal,NAdtotal,18,18) :: gout
  COMPLEX(double), DIMENSION(Nadtotal,Nadtotal,18,18) :: gtemp

  COMPLEX(double), DIMENSION(2,Nadtotal,Nadtotal,9,9) :: g_upup,g_updn,g_dnup,g_dndn

  COMPLEX(double) :: wc

  INTEGER :: a,b,ll,l,m,n,mu,mul,nu,nul

  wc = CMPLX(csi+w,eta)
  CALL green(wc,gtemp)
  gout(1,:,:,:,:) = gtemp

  wc = CMPLX(csi,eta)
  CALL green(wc,gtemp)
  gout(2,:,:,:,:) = gtemp

  DO l=1,2
  DO is=1,Nadtotal
  DO js=1,Nadtotal
  g_upup(l,is,js,:,:) = gout(l,is,js,1:9,1:9)
  g_updn(l,is,js,:,:) = gout(l,is,js,1:9,10:18)
  g_dnup(l,is,js,:,:) = gout(l,is,js,10:18,1:9)
  g_dndn(l,is,js,:,:) = gout(l,is,js,10:18,10:18)
  END DO
  END DO
  END DO

  l = 1
  DO is=1,Nadtotal
  DO js=1,Nadtotal
  DO mu = 1,5
  DO mul = 1,5
  m = (l-1)*25*Nadtotal +  (is-1)*25 + (mu-1)*5 + mul 
  DO nu = 1,5
  DO nul = 1,5
  n = (js-1)*25 + (nu-1)*5 + nul

  chi0m(m,n) =  g_dndn(1,is,js,nu+4,mul+4)*CONJG(g_upup(2,is,js,mu+4,nul+4))
  END DO
  END DO
  END DO
  END DO
  END DO
  END DO

  l = 2
  DO is=1,Nadtotal
  DO js=1,Nadtotal
  DO mu = 1,5
  DO mul = 1,5
  m = (l-1)*25*Nadtotal + (is-1)*25 + (mu-1)*5 + mul 
  DO nu = 1,5
  DO nul = 1,5
  n = (js-1)*25 + (nu-1)*5 + nul

  chi0m(m,n) =  g_updn(1,is,js,nu+4,mul+4)*CONJG(g_upup(2,is,js,mu+4,nul+4))
  END DO
  END DO
  END DO
  END DO
  END DO
  END DO

  l = 3
  DO is=1,Nadtotal
  DO js=1,Nadtotal
  DO mu = 1,5
  DO mul = 1,5
  m = (l-1)*25*Nadtotal + (is-1)*25 + (mu-1)*5 + mul 
  DO nu = 1,5
  DO nul = 1,5
  n = (js-1)*25 + (nu-1)*5 + nul

  chi0m(m,n) =  g_dndn(1,is,js,nu+4,mul+4)*CONJG(g_dnup(2,is,js,mu+4,nul+4))
  END DO
  END DO
  END DO
  END DO
  END DO
  END DO
 
  l = 4
  DO is=1,Nadtotal
  DO js=1,Nadtotal
  DO mu = 1,5
  DO mul = 1,5
  m = (l-1)*25*Nadtotal + (is-1)*25 + (mu-1)*5 + mul 
  DO nu = 1,5
  DO nul = 1,5
  n = (js-1)*25 + (nu-1)*5 + nul

  chi0m(m,n) = g_updn(1,is,js,nu+4,mul+4)*CONJG(g_dnup(2,is,js,mu+4,nul+4))
  END DO
  END DO
  END DO
  END DO
  END DO
  END DO

  chi0m = -zi*chi0m/(2.d0*pi)

  RETURN

END SUBROUTINE dchi_HF_2





SUBROUTINE dLambda_1(w,csi,Lambdam)

 USE f90_kind
 USE constants_u
 USE ef_and_eta
 USE lattice
 USE MPI_pars


  REAL(double), INTENT(IN) :: w,csi
  COMPLEX(double), DIMENSION(4*5*5*Nadtotal,4*5*5*Nadtotal), INTENT(OUT) :: Lambdam

  COMPLEX(double), DIMENSION(4*5*5*Nadtotal,4*5*5*Nadtotal) :: Lambda,dLambda

  COMPLEX(double), DIMENSION(4,4,Nadtotal,Nadtotal,5,5) :: Kappa

  COMPLEX(double), DIMENSION(3,Nadtotal,Nadtotal,18,18) :: gout
  COMPLEX(double), DIMENSION(Nadtotal,Nadtotal,18,18) :: gtemp

  COMPLEX(double), DIMENSION(3,Nadtotal,Nadtotal,9,9) :: g_upup,g_updn,g_dnup,g_dndn

  COMPLEX(double) :: wc

  INTEGER :: a,b,a0,b0,l,ll,m,n,mu,mul,nu,nul


  Lambda = zero
  dLambda = zero

  wc = CMPLX(ef+w,csi)
  CALL green(wc,gtemp)
  gout(1,:,:,:,:) = gtemp

  wc = CMPLX(ef,csi)
  CALL green(wc,gtemp)
  gout(2,:,:,:,:) = gtemp
  
  wc = CMPLX(ef-w,csi)
  CALL green(wc,gtemp)
  gout(3,:,:,:,:) = gtemp


  DO l=1,3
  DO is=1,Nadtotal
  DO js=1,Nadtotal
  g_upup(l,is,js,:,:) = gout(l,is,js,1:9,1:9)
  g_updn(l,is,js,:,:) = gout(l,is,js,1:9,10:18)
  g_dnup(l,is,js,:,:) = gout(l,is,js,10:18,1:9)
  g_dndn(l,is,js,:,:) = gout(l,is,js,10:18,10:18)
  END DO
  END DO
  END DO

  Kappa = zero
  DO is=1,Nadtotal
  DO js=1,Nadtotal
  DO mu = 5,9
  DO nu = 5,9
  DO nul = 5,9
  Kappa(1,1,is,js,mu-4,nu-4) = Kappa(1,1,is,js,mu-4,nu-4) - &
                               g_dndn(1,is,js,nu,nul)*g_upup(2,js,is,nul,mu) - &
                               CONJG(g_upup(3,is,js,mu,nul)*g_dndn(2,js,is,nul,nu))

  Kappa(1,2,is,js,mu-4,nu-4) = Kappa(1,2,is,js,mu-4,nu-4) - & 
                               g_dnup(1,is,js,nu,nul)*g_upup(2,js,is,nul,mu) -       &
                               CONJG(g_upup(3,is,js,mu,nul)*g_updn(2,js,is,nul,nu)) 

  Kappa(1,3,is,js,mu-4,nu-4) = Kappa(1,3,is,js,mu-4,nu-4) -              &
                               g_dndn(1,is,js,nu,nul)*g_dnup(2,js,is,nul,mu) - &
                               CONJG(g_updn(3,is,js,mu,nul)*g_dndn(2,js,is,nul,nu))

  Kappa(1,4,is,js,mu-4,nu-4) = Kappa(1,4,is,js,mu-4,nu-4) -              &
                               g_dnup(1,is,js,nu,nul)*g_dnup(2,js,is,nul,mu) - &
                               CONJG(g_updn(3,is,js,mu,nul)*g_updn(2,js,is,nul,nu))
 
  Kappa(2,1,is,js,mu-4,nu-4) = Kappa(2,1,is,js,mu-4,nu-4) -              &
                               g_updn(1,is,js,nu,nul)*g_upup(2,js,is,nul,mu) - &
                               CONJG(g_upup(3,is,js,mu,nul)*g_dnup(2,js,is,nul,nu))

  Kappa(2,2,is,js,mu-4,nu-4) = Kappa(2,2,is,js,mu-4,nu-4) -              &
                               CONJG(g_upup(3,is,js,mu,nul)*g_upup(2,js,is,nul,nu)) - &
                               g_upup(1,is,js,nu,nul)*g_upup(2,js,is,nul,mu)

  Kappa(2,3,is,js,mu-4,nu-4) = Kappa(2,3,is,js,mu-4,nu-4) - &
                               g_updn(1,is,js,nu,nul)*g_dnup(2,js,is,nul,mu) - &
                               CONJG(g_updn(3,is,js,mu,nul)*g_dnup(2,js,is,nul,nu))

  Kappa(2,4,is,js,mu-4,nu-4) = Kappa(2,4,is,js,mu-4,nu-4) - &
                               g_upup(1,is,js,nu,nul)*g_dnup(2,js,is,nul,mu) - &
                               CONJG(g_updn(3,is,js,mu,nul)*g_upup(2,js,is,nul,nu))

  Kappa(3,1,is,js,mu-4,nu-4) = Kappa(3,1,is,js,mu-4,nu-4) - &
                               g_dndn(1,is,js,nu,nul)*g_updn(2,js,is,nul,mu) - &
                               CONJG(g_dnup(3,is,js,mu,nul)*g_dndn(2,js,is,nul,nu))

  Kappa(3,2,is,js,mu-4,nu-4) = Kappa(3,2,is,js,mu-4,nu-4) - &
                               g_dnup(1,is,js,nu,nul)*g_updn(2,js,is,nul,mu) - &
                               CONJG(g_dnup(3,is,js,mu,nul)*g_updn(2,js,is,nul,nu))

  Kappa(3,3,is,js,mu-4,nu-4) = Kappa(3,3,is,js,mu-4,nu-4) - &
                               g_dndn(1,is,js,nu,nul)*g_dndn(2,js,is,nul,mu) - &
                               CONJG(g_dndn(3,is,js,mu,nul)*g_dndn(2,js,is,nul,nu))

  Kappa(3,4,is,js,mu-4,nu-4) = Kappa(3,4,is,js,mu-4,nu-4) - &
                               g_dnup(1,is,js,nu,nul)*g_dndn(2,js,is,nul,mu) - &
                               CONJG(g_dndn(3,is,js,mu,nul)*g_updn(2,js,is,nul,nu))

  Kappa(4,1,is,js,mu-4,nu-4) = Kappa(4,1,is,js,mu-4,nu-4) - &
                               g_updn(1,is,js,nu,nul)*g_updn(2,js,is,nul,mu) - &
                               CONJG(g_dnup(3,is,js,mu,nul)*g_dnup(2,js,is,nul,nu))

  Kappa(4,2,is,js,mu-4,nu-4) = Kappa(4,2,is,js,mu-4,nu-4) - &
                               g_upup(1,is,js,nu,nul)*g_updn(2,js,is,nul,mu) - &
                               CONJG(g_dnup(3,is,js,mu,nul)*g_upup(2,js,is,nul,nu))

  Kappa(4,3,is,js,mu-4,nu-4) = Kappa(4,3,is,js,mu-4,nu-4) - &
                               g_updn(1,is,js,nu,nul)*g_dndn(2,js,is,nul,mu) - &
                               CONJG(g_dndn(3,is,js,mu,nul)*g_dnup(2,js,is,nul,nu))

  Kappa(4,4,is,js,mu-4,nu-4) = Kappa(4,4,is,js,mu-4,nu-4) - &
                               g_upup(1,is,js,nu,nul)*g_dndn(2,js,is,nul,mu) - &
                               CONJG(g_dndn(3,is,js,mu,nul)*g_upup(2,js,is,nul,nu))
  END DO
  END DO
  END DO
  END DO
  END DO

      
  DO is=1,Nadtotal
  DO js=1,Nadtotal
  DO mu = 5,9
  DO nu = 5,9
  DO mul = 5,9

  DO nul = 5,9

  l = 1
  ll = 2
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + nul - 4
  a0 = a
  b0 = b
  dLambda(a,b) = g_dndn(1,is,js,nu,nul)*g_dnup(2,js,is,mul,mu) + &
                 CONJG(g_upup(3,is,js,mu,mul)*g_updn(2,js,is,nul,nu)) + &
                 g_dnup(1,is,js,nu,nul)*g_upup(2,js,is,mul,mu) +        &
                 CONJG(g_updn(3,is,js,mu,mul)*g_dndn(2,js,is,nul,nu))

  l = 1
  ll = 3
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + nul - 4
  dLambda(a,b) = dLambda(a0,b0)

  l = 2
  ll = 2
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + nul - 4
  a0 = a
  b0 = b
  dLambda(a,b) = g_updn(1,is,js,nu,nul)*g_dnup(2,js,is,mul,mu) + &
                 CONJG(g_upup(3,is,js,mu,mul)*g_upup(2,js,is,nul,nu)) + &
                 g_upup(1,is,js,nu,nul)*g_upup(2,js,is,mul,mu) +        &
                 CONJG(g_updn(3,is,js,mu,mul)*g_dnup(2,js,is,nul,nu))
  
  l = 2
  ll = 3
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + nul - 4
  dLambda(a,b)  = dLambda(a0,b0) ! = dLambda(2,2) except for Kappa(2,3)
   
  l = 3
  ll = 2
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + nul - 4
  a0 = a
  b0 = b
  dLambda(a,b) = g_dndn(1,is,js,nu,nul)*g_dndn(2,js,is,mul,mu) + &
                 CONJG(g_dnup(3,is,js,mu,mul)*g_updn(2,js,is,nul,nu)) + &
                 g_dnup(1,is,js,nu,nul)*g_updn(2,js,is,mul,mu) + &
                 CONJG(g_dndn(3,is,js,mu,mul)*g_dndn(2,js,is,nul,nu))
  
  l = 3
  ll = 3
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + nul - 4
  dLambda(a,b)  = dLambda(a0,b0) ! = dLambda(3,2) except for Kappa(3,3)
   
  l = 4
  ll = 2
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + nul - 4
  a0 = a
  b0 = b
  dLambda(a,b) = g_updn(1,is,js,nu,nul)*g_dndn(2,js,is,mul,mu) + &
                 CONJG(g_dnup(3,is,js,mu,mul)*g_upup(2,js,is,nul,nu)) + &
                 g_upup(1,is,js,nu,nul)*g_updn(2,js,is,mul,mu) + &
                 CONJG(g_dndn(3,is,js,mu,mul)*g_dnup(2,js,is,nul,nu))
  
  l = 4
  ll = 3
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + nul - 4
  dLambda(a,b) = dLambda(a0,b0) ! = dLambda(4,2) except for Kappa(3,3)
  
  END DO ! nul=5,9

  l = 1
  ll = 1
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = Kappa(1,1,is,js,mu-4,nu-4)
 
  l = 1
  ll = 2
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = dLambda(a,b) + Kappa(1,2,is,js,mu-4,nu-4)
 
  l = 1
  ll = 3
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = dLambda(a,b) + Kappa(1,3,is,js,mu-4,nu-4)
 
  l = 1
  ll = 4
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = Kappa(1,4,is,js,mu-4,nu-4)
 
  l = 2
  ll = 1
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = Kappa(2,1,is,js,mu-4,nu-4)
 
  l = 2
  ll = 2
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = dLambda(a,b) + Kappa(2,2,is,js,mu-4,nu-4)
 
  l = 2
  ll = 3
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = dLambda(a,b) + Kappa(2,3,is,js,mu-4,nu-4)

  l = 2
  ll = 4
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = Kappa(2,4,is,js,mu-4,nu-4)
  
  l = 3
  ll = 1
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = Kappa(3,1,is,js,mu-4,nu-4)

  l = 3
  ll = 2
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = dLambda(a,b) + Kappa(3,2,is,js,mu-4,nu-4)
  
  l = 3
  ll = 3
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = dLambda(a,b) + Kappa(3,3,is,js,mu-4,nu-4)
   
  l = 3
  ll = 4
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = Kappa(3,4,is,js,mu-4,nu-4)
   
  l = 4
  ll = 1
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = Kappa(4,1,is,js,mu-4,nu-4)
 
  l = 4
  ll = 2
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = dLambda(a,b) + Kappa(4,2,is,js,mu-4,nu-4)
 
  l = 4
  ll = 3
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = dLambda(a,b) + Kappa(4,3,is,js,mu-4,nu-4)

  l = 4
  ll = 4
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = Kappa(4,4,is,js,mu-4,nu-4)


  END DO ! mul=5,9
  END DO ! nu=5,9
  END DO ! mu=5,9
  END DO ! js=1,Nadtotal
  END DO ! is=1,Nadtotal

  Lambdam = dLambda/(2.d0*pi)


  RETURN

END SUBROUTINE dLambda_1





SUBROUTINE dLambda_2(w,csi,Lambdam)

 USE f90_kind
 USE constants_u
 USE ef_and_eta
 USE lattice
 USE MPI_pars

  REAL(double), INTENT(IN) :: w,csi

  COMPLEX(double), DIMENSION(4*5*5*Nadtotal,4*5*5*Nadtotal), INTENT(OUT) :: Lambdam

  COMPLEX(double), DIMENSION(4*5*5*Nadtotal,4*5*5*Nadtotal) :: dLambda

  COMPLEX(double), DIMENSION(4,4,Nadtotal,Nadtotal,5,5) :: Kappa

  COMPLEX(double), DIMENSION(2,Nadtotal,Nadtotal,18,18) :: gout
  COMPLEX(double), DIMENSION(Nadtotal,Nadtotal,18,18) :: gtemp

  COMPLEX(double), DIMENSION(2,Nadtotal,Nadtotal,9,9) :: g_upup,g_updn,g_dnup,g_dndn

  COMPLEX(double) :: wc

  INTEGER :: a,b,a0,b0,ll,l,m,n,mu,mul,nu,nul


  dLambda = zero

  wc = CMPLX(csi+w,eta)
  CALL green(wc,gtemp)
  gout(1,:,:,:,:) = gtemp

  wc = CMPLX(csi,eta)
  CALL green(wc,gtemp)
  gout(2,:,:,:,:) = gtemp
  
  DO l=1,2
  DO is=1,Nadtotal
  DO js=1,Nadtotal
  g_upup(l,is,js,:,:) = gout(l,is,js,1:9,1:9)
  g_updn(l,is,js,:,:) = gout(l,is,js,1:9,10:18)
  g_dnup(l,is,js,:,:) = gout(l,is,js,10:18,1:9)
  g_dndn(l,is,js,:,:) = gout(l,is,js,10:18,10:18)
  END DO
  END DO
  END DO

  Kappa = zero
  DO is=1,Nadtotal
  DO js=1,Nadtotal
  DO mu = 5,9
  DO nu = 5,9
  DO nul = 5,9
  Kappa(1,1,is,js,mu-4,nu-4) = Kappa(1,1,is,js,mu-4,nu-4) - &
  g_dndn(1,is,js,nu,nul)*CONJG(g_upup(2,is,js,mu,nul))

  Kappa(1,2,is,js,mu-4,nu-4) = Kappa(1,2,is,js,mu-4,nu-4) -                 & 
  g_dnup(1,is,js,nu,nul)*CONJG(g_upup(2,is,js,mu,nul))

  Kappa(1,3,is,js,mu-4,nu-4) = Kappa(1,3,is,js,mu-4,nu-4) -              &
  g_dndn(1,is,js,nu,nul)*CONJG(g_updn(2,is,js,mu,nul))

  Kappa(1,4,is,js,mu-4,nu-4) = Kappa(1,4,is,js,mu-4,nu-4) -              &
  g_dnup(1,is,js,nu,nul)*CONJG(g_updn(2,is,js,mu,nul))

  Kappa(2,1,is,js,mu-4,nu-4) = Kappa(2,1,is,js,mu-4,nu-4) -              &
  g_updn(1,is,js,nu,nul)*CONJG(g_upup(2,is,js,mu,nul))

  Kappa(2,2,is,js,mu-4,nu-4) = Kappa(2,2,is,js,mu-4,nu-4) -              &
  g_upup(1,is,js,nu,nul)*CONJG(g_upup(2,is,js,mu,nul))

  Kappa(2,3,is,js,mu-4,nu-4) = Kappa(2,3,is,js,mu-4,nu-4) - &
  g_updn(1,is,js,nu,nul)*CONJG(g_updn(2,is,js,mu,nul))

  Kappa(2,4,is,js,mu-4,nu-4) = Kappa(2,4,is,js,mu-4,nu-4) - &
  g_upup(1,is,js,nu,nul)*CONJG(g_updn(2,is,js,mu,nul))

  Kappa(3,1,is,js,mu-4,nu-4) = Kappa(3,1,is,js,mu-4,nu-4) - &
  g_dndn(1,is,js,nu,nul)*CONJG(g_dnup(2,is,js,mu,nul))

  Kappa(3,2,is,js,mu-4,nu-4) = Kappa(3,2,is,js,mu-4,nu-4) - &
  g_dnup(1,is,js,nu,nul)*CONJG(g_dnup(2,is,js,mu,nul)) 

  Kappa(3,3,is,js,mu-4,nu-4) = Kappa(3,3,is,js,mu-4,nu-4) - &
  g_dndn(1,is,js,nu,nul)*CONJG(g_dndn(2,is,js,mu,nul)) 

  Kappa(3,4,is,js,mu-4,nu-4) = Kappa(3,4,is,js,mu-4,nu-4) - &
  g_dnup(1,is,js,nu,nul)*CONJG(g_dndn(2,is,js,mu,nul)) 

  Kappa(4,1,is,js,mu-4,nu-4) = Kappa(4,1,is,js,mu-4,nu-4) - &
  g_updn(1,is,js,nu,nul)*CONJG(g_dnup(2,is,js,mu,nul)) 

  Kappa(4,2,is,js,mu-4,nu-4) = Kappa(4,2,is,js,mu-4,nu-4) - &
  g_upup(1,is,js,nu,nul)*CONJG(g_dnup(2,is,js,mu,nul))

  Kappa(4,3,is,js,mu-4,nu-4) = Kappa(4,3,is,js,mu-4,nu-4) - &
  g_updn(1,is,js,nu,nul)*CONJG(g_dndn(2,is,js,mu,nul)) 

  Kappa(4,4,is,js,mu-4,nu-4) = Kappa(4,4,is,js,mu-4,nu-4) - &
  g_upup(1,is,js,nu,nul)*CONJG(g_dndn(2,is,js,mu,nul)) 

  END DO
  END DO
  END DO
  END DO
  END DO


  DO is=1,Nadtotal
  DO js=1,Nadtotal
  DO mu = 5,9
  DO nu = 5,9
  DO mul = 5,9

  DO nul = 5,9

  l = 1
  ll= 2
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + nul - 4
  a0 = a
  b0 = b
  dLambda(a,b) = g_dndn(1,is,js,nu,nul)*CONJG(g_updn(2,is,js,mu,mul)) + &
                 g_dnup(1,is,js,nu,nul)*CONJG(g_upup(2,is,js,mu,mul))

  l = 1
  ll= 3
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + nul - 4
  dLambda(a,b) = dLambda(a0,b0)! = dLambda(1,2) ** (except for Kappa)


  l = 2
  ll= 2
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + nul - 4
  a0 = a
  b0 = b
  dLambda(a,b) = g_updn(1,is,js,nu,nul)*CONJG(g_updn(2,is,js,mu,mul)) + &
                 g_upup(1,is,js,nu,nul)*CONJG(g_upup(2,is,js,mu,mul)) 

  l = 2
  ll= 3
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + nul - 4
  dLambda(a,b) = dLambda(a0,b0) ! = dLambda(2,2) except for Kappa(2,3)

  l = 3
  ll= 2
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + nul - 4
  a0 = a
  b0 = b
  dLambda(a,b) = g_dndn(1,is,js,nu,nul)*CONJG(g_dndn(2,is,js,mu,mul)) + &
                 g_dnup(1,is,js,nu,nul)*CONJG(g_dnup(2,is,js,mu,mul)) 

  l = 3
  ll= 3
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + nul - 4
  dLambda(a,b) = dLambda(a0,b0) ! = dLambda(3,2) except for Kappa(3,3)

  l = 4
  ll= 2
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + nul - 4
  a0 = a
  b0 = b
  dLambda(a,b) = g_updn(1,is,js,nu,nul)*CONJG(g_dndn(2,is,js,mu,mul)) + &
                 g_upup(1,is,js,nu,nul)*CONJG(g_dnup(2,is,js,mu,mul))

  l = 4
  ll= 3
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + nul - 4
  dLambda(a,b) = dLambda(a0,b0)! = dLambda(4,2) except for Kappa(3,3)

  END DO ! nul=5,9

  l = 1
  ll= 1
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = Kappa(1,1,is,js,mu-4,nu-4)

  l = 1
  ll= 2
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = dLambda(a,b) + Kappa(1,2,is,js,mu-4,nu-4)
                       
  l = 1
  ll= 3
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = dLambda(a,b) + Kappa(1,3,is,js,mu-4,nu-4)
                       
  l = 1
  ll= 4
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = Kappa(1,4,is,js,mu-4,nu-4)
                       
  l = 2
  ll= 1
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = Kappa(2,1,is,js,mu-4,nu-4)
                       
  l = 2
  ll= 2
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = dLambda(a,b) + Kappa(2,2,is,js,mu-4,nu-4)
                       
  l = 2
  ll= 3
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = dLambda(a,b) + Kappa(2,3,is,js,mu-4,nu-4)
                       
  l = 2
  ll= 4
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = Kappa(2,4,is,js,mu-4,nu-4)
                       
  l = 3
  ll= 1
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = Kappa(3,1,is,js,mu-4,nu-4)
                       
  l = 3
  ll= 2
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = dLambda(a,b) + Kappa(3,2,is,js,mu-4,nu-4)
                       
  l = 3
  ll= 3
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = dLambda(a,b) + Kappa(3,3,is,js,mu-4,nu-4)
                       
  l = 3
  ll= 4
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = Kappa(3,4,is,js,mu-4,nu-4)
                       
  l = 4
  ll= 1
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = Kappa(4,1,is,js,mu-4,nu-4)
                       
  l = 4
  ll= 2
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = dLambda(a,b) + Kappa(4,2,is,js,mu-4,nu-4)
                       
  l = 4
  ll= 3
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = dLambda(a,b) + Kappa(4,3,is,js,mu-4,nu-4)
                       
  l = 4
  ll= 4
  a = (l-1)*25*Nadtotal + (is-1)*25 + 5*(mu-5) + nu - 4
  b = (ll-1)*25*Nadtotal + (js-1)*25 + 5*(mul-5) + mul - 4
  dLambda(a,b) = Kappa(4,4,is,js,mu-4,nu-4)


  END DO ! mul=5,9
  END DO ! nu=5,9
  END DO ! mu=5,9
  END DO ! js=1,Nadtotal
  END DO ! is=1,Nadtotal

  Lambdam = -zi*dLambda/(2.d0*pi)


  RETURN

END SUBROUTINE dLambda_2
