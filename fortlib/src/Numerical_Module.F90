module Numerical_Module

  implicit none

contains
! compute eigenvalue & eigenvector by Jabobi method
subroutine Jacobi_method(A,X,n)
! state variables
implicit none
integer :: i,k,l
integer, intent(in) :: n
integer :: p,q,fg
real(4) :: bg,ss,tt,vv,sai,sg,csai
real(4) :: bpi,bqi,bpp,bqq,xpq,xqp
real(4) :: eps
real(4), intent(in) :: A(n,n)
real(4), intent(out) :: X(n,n)
! n is the degree of matrix (You should check.)
parameter (n=5)
allocate(A(n,n))
allocate(X(n,n))

do l=1,n
 write(*,'(16(E10.3,X))') (A(l,k),k=1,n)
end do
! judge if data is symmetric matrix
do l=1,n
 do k=1,n
  if (A(l,k).EQ.A(k,l)) then
  else
   write(*,200)
   200 format('It is not symmetric matrix.')
  end if
 end do
end do
! set initial value of matrix X
do l=1,n
 do k=1,n
  X(l,k)=0
 end do
 X(l,l)=1
end do
! find maximum element of matrix A
fg=0
bg=0
do k=1,n-1
 do l=k+1,n
  if (abs(a(k,l)).GT.bg) then
   p=k
   q=l
   bg=abs(a(k,l))
  end if
 end do
end do
! compute orthogonal matrix & diagonalize matrix A
do i=1,1000
 if (bg.GE.eps) then
  fg=fg+1
  ss=-A(p,q)
  tt=(A(p,p)-A(q,q))/2
  vv=abs(tt)/sqrt(ss*ss+tt*tt)
  sai=sqrt((1-vv)/2)
  sg=ss*tt
  csai=sqrt(1-sai*sai)
  if (sg.LT.0) then
   sai=-sai
  end if
  do l=1,n
   if (l.NE.p.AND.l.NE.q) then
    bpi=A(p,l)*csai-A(q,l)*sai
    bqi=A(q,l)*csai+A(p,l)*sai
    A(p,l)=bpi
    A(q,l)=bqi
   end if
  end do
  bpp=A(p,p)*csai*csai+A(q,q)*sai*sai-2*A(p,q)*sai*csai
  bqq=A(p,p)*sai*sai+A(q,q)*csai*csai+2*A(p,q)*sai*csai
  A(p,p)=bpp
  A(q,q)=bqq
  A(p,q)=0
  A(q,p)=0
! compute eigenvector
  do l=1,n
   A(l,p)=A(p,l)
   A(l,q)=A(q,l)
   xpq=X(l,p)*csai-X(l,q)*sai
   xqp=X(l,q)*csai+X(l,p)*sai
   X(l,p)=xpq
   X(l,q)=xqp
  end do
! find maximum element of matrix A
  bg=0
  do k=1,n-1
   do l=k+1,n
    if (abs(A(k,l)).GT.bg) then
     p=k
     q=l
     bg=abs(A(k,l))
    end if
   end do
  end do
 end if
end do

end subroutine Jacobi_method
