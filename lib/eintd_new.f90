subroutine eint_spd(sk,w,b)

 use f90_kind

 REAL(double), DIMENSION(10), INTENT(IN) :: sk
 REAL(double), DIMENSION(3), INTENT(IN) :: w
 REAL(double), DIMENSION(9,9), INTENT(OUT) :: b

 REAL(double) :: sss,sps,pps,ppp,ss,ps,pp,ds,dp,dd
 REAL(double) :: x,y,z,xx,xy,yy,yz,zz,zx,xxyy,yyzz,zzxx,aux,r3,      &
                 aux1,f8,f1,f2,f3,g1,g2,f4,f5,aux2,aux3,aux4   


  sss = sk(1)
  sps = sk(2)
  pps = sk(3)
  ppp = sk(4)
  ss = sk(5)
  ps = sk(6)
  pp = sk(7)
  ds = sk(8)
  dp = sk(9)
  dd = sk(10)

  x=w(1)
  y=w(2)
  z=w(3)

  xx=x*x
  xy=x*y
  yy=y*y
  yz=y*z
  zz=z*z
  zx=z*x
  xxyy=xx*yy
  yyzz=yy*zz
  zzxx=zz*xx
  aux=pps-ppp
  r3=dsqrt(3.d0)
  aux1=r3*ss
  f8=3.d0*zz-1.d0
  f1=xx+yy
  f2=xx-yy
  f3=zz-.5d0*f1
  g1=1.5d0*f2*ds
  g2=r3*f3*ds
  b(1,1)=sss
  b(1,2)=x*sps
  b(1,3)=y*sps
  b(1,4)=z*sps
  b(2,1)=-b(1,2)
  b(2,2)=xx*pps+(1.d0-xx)*ppp
  b(2,3)=xy*aux
  b(2,4)=zx*aux
  b(3,1)=-b(1,3)
  b(3,2)=b(2,3)
  b(3,3)=yy*pps+(1.d0-yy)*ppp
  b(3,4)=yz*aux
  b(4,1)=-b(1,4)
  b(4,2)=b(2,4)
  b(4,3)=b(3,4)
  b(4,4)=zz*pps+(1.d0-zz)*ppp
  b(1,5)=xy*aux1
  b(1,6)=yz*aux1
  b(1,7)=zx*aux1
  b(1,8)=.5d0*f2*aux1
  b(1,9)=.5d0*f8*ss
  b(5,1)=b(1,5)
  b(6,1)=b(1,6)
  b(7,1)=b(1,7)
  b(8,1)=b(1,8)
  b(9,1)=b(1,9)
  f4=.5d0*r3*f2*ps
  f5=.5d0*f8*ps
  aux2=r3*xx*ps+(1.d0-2.d0*xx)*pp
  b(2,5)=aux2*y
  b(2,6)=(r3*ps-2.d0*pp)*xy*z
  b(2,7)=aux2*z
  b(2,8)=(f4+(1.d0-f2)*pp)*x
  b(2,9)=(f5-r3*zz*pp)*x
  aux3=(r3*yy*ps+(1.d0-2.d0*yy)*pp)
  b(3,5)=aux3*x
  b(3,6)=aux3*z
  b(3,7)=b(2,6)
  b(3,8)=(f4-(1.d0+f2)*pp)*y
  b(3,9)=(f5-r3*zz*pp)*y
  aux4=r3*zz*ps+(1.d0-2.d0*zz)*pp
  b(4,5)=b(2,6)
  b(4,6)=aux4*y
  b(4,7)=aux4*x
  b(4,8)=(f4-f2*pp)*z
  b(4,9)=(f5+r3*f1*pp)*z
  b(5,2)=-b(2,5)
  b(6,2)=-b(2,6)
  b(7,2)=-b(2,7)
  b(8,2)=-b(2,8)
  b(9,2)=-b(2,9)
  b(5,3)=-b(3,5)
  b(6,3)=-b(3,6)
  b(7,3)=-b(3,7)
  b(8,3)=-b(3,8)
  b(9,3)=-b(3,9)
  b(5,4)=-b(4,5)
  b(6,4)=-b(4,6)
  b(7,4)=-b(4,7)
  b(8,4)=-b(4,8)
  b(9,4)=-b(4,9)
  b(5,5)=3.d0*xxyy*ds+(f1-4.d0*xxyy)*dp+(zz+xxyy)*dd
  b(5,6)=(3.d0*yy*ds+(1.d0-4.d0*yy)*dp+(yy-1.d0)*dd)*zx
  b(5,7)=(3.d0*xx*ds+(1.d0-4.d0*xx)*dp+(xx-1.d0)*dd)*yz
  b(5,8)=(g1-2.d0*f2*dp+.5d0*f2*dd)*xy
  b(5,9)=(g2-2.d0*r3*zz*dp+.5d0*r3*(1.d0+zz)*dd)*xy
  b(6,5)=b(5,6)
  b(6,6)=3.d0*yyzz*ds+(yy+zz-4.d0*yyzz)*dp+(xx+yyzz)*dd
  b(6,7)=(3.d0*zz*ds+(1.d0-4.d0*zz)*dp+(zz-1.d0)*dd)*xy
  b(6,8)=(g1-(1.d0+2.d0*f2)*dp+(1.d0+.5d0*f2)*dd)*yz
  b(6,9)=(g2+r3*(f1-zz)*dp-.5d0*r3*f1*dd)*yz
  b(7,5)=b(5,7)
  b(7,6)=b(6,7)
  b(7,7)=3.d0*zzxx*ds+(zz+xx-4.d0*zzxx)*dp+(yy+zzxx)*dd
  b(7,8)=(g1+(1.d0-2.d0*f2)*dp-(1.d0-.5d0*f2)*dd)*zx
  b(7,9)=(g2+r3*(f1-zz)*dp-.5d0*r3*f1*dd)*zx
  b(8,5)=b(5,8)
  b(8,6)=b(6,8)
  b(8,7)=b(7,8)
  b(8,8)=.75d0*f2*f2*ds+(f1-f2*f2)*dp+(zz+.25d0*f2*f2)*dd
  b(8,9)=.5d0*f2*g2-r3*zz*f2*dp+.25d0*r3*(1.d0+zz)*f2*dd
  b(9,5)=b(5,9)
  b(9,6)=b(6,9)
  b(9,7)=b(7,9)
  b(9,8)=b(8,9)
  b(9,9)=f3*f3*ds+3.d0*zz*f1*dp+.75d0*f1*f1*dd

  return

end subroutine eint_spd
