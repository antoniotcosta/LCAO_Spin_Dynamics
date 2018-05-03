subroutine eint(sk,w,b)

 use f90_kind

 REAL(double), DIMENSION(4), INTENT(IN) :: sk
 REAL(double), DIMENSION(3), INTENT(IN) :: w
 REAL(double), DIMENSION(4,4), INTENT(OUT) :: b

 REAL(double) :: sss,sps,pps,ppp
 REAL(double) :: x,y,z,xx,xy,yy,yz,zz,zx,xxyy,yyzz,zzxx,aux,r3,      &
                 aux1,f8,f1,f2,f3,g1,g2,f4,f5,aux2,aux3,aux4   


  sss = sk(1)
  sps = sk(2)
  pps = sk(3)
  ppp = sk(4)

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
  f8=3.d0*zz-1.d0
  f1=xx+yy
  f2=xx-yy
  f3=zz-.5d0*f1
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

  return

end subroutine eint
