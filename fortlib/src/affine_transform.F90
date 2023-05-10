module AFFINE_TRANSFORMATION
  contains
  

  SUBROUTINE CALC_AFFINE_COEFF_FOR_4_POINTS(before,after,af)
!    USE f95_lapack
    USE f77_lapack
    IMPLICIT NONE

  integer :: i, j,ii,jj,info, small_mn, nmmax,inum
  integer,parameter :: m=8
  integer,parameter :: n=6
  REAL(8) :: x1,x2,y1,y2,xx,yy,xy  
  REAL(4),INTENT(IN) :: before(2,4),after(2,4)
  REAL(4),INTENT(OUT) :: Af(n)
  real(8) :: A(m,n)
  real(8) :: B(m)
  real(8),allocatable :: S(:), Si(:), U(:,:), Vt(:,:),Ut(:,:)
  real(8),allocatable :: dum1(:,:),dum2(:,:),dum3(:,:)
  real(8) :: Sd 
  integer,parameter :: lwork = 50
  real(8) :: work(50)
  
  inum=0
  do j =1,4
     do i =1,2
        inum=1+inum
     b(inum)=after(i,j)
  enddo
enddo

a(:,:)=0
  

inum=0
  
do j =1,m
   if(mod(j,2).ne.0) then
    jj=(j+1)/2;a(j,1:2)=before(1:2,jj);a(j,3)=1
    else
    jj=j/2;a(j,4:5)=before(1:2,jj);a(j,6)=1
    endif
 enddo

!!$ do i =1,m
!!$    print *,b(i),a(i,:)
!!$ enddo
 
  small_mn = min(m,n)
  allocate(S(1:small_mn),Si(1:small_mn), U(1:m,1:m),Vt(1:n,1:n))
  allocate(dum1(n,m),dum3(n,m),dum2(n,m),Ut(1:m,1:m))
  
  call dgesvd('A','A',m,n,A(1:m,1:n), m, S(1:small_mn), U,m, &
        Vt,n,work,lwork,info )
!  call LA_GESVD(A(1:m,1:n), S(1:small_mn), U=U(1:m,1:m), &
!        Vt = Vt(1:n,1:n),JOB ="N", INFO = info )

  if(info < 0) then 
     write(*,*) info
     stop "Fail"
  endif
  
  Sd=product(S)
  Si(:)=1
  do j=1,small_mn
     do i=1,small_mn
        if(j.eq.i) cycle
        Si(j)=S(i)*Si(j)
     enddo
  enddo
  
  Si(:)=Si(:)/Sd
  
  dum3(:,:)=0
  
  do i =1,small_mn
     dum3(i,i)=Si(i)
  enddo
  
  Ut=transpose(U)
  Vt=transpose(Vt)  
  dum1=matmul(VT,dum3)
  dum2=matmul(dum1,Ut)
  af=matmul(dum2,B)
  
END SUBROUTINE CALC_AFFINE_COEFF_FOR_4_POINTS


subroutine TRANSFORMATION(x,y,af,x2,y2)
  IMPLICIT NONE
  REAL(4),INTENT(IN) :: x,y,af(6)
  REAL(4),INTENT(OUT) :: x2,y2

  x2=af(1)*x+af(2)*y+af(3)
  y2=af(4)*x+af(5)*y+af(6)
endsubroutine TRANSFORMATION
  

endmodule AFFINE_TRANSFORMATION
