module Numerical_Module
  IMPLICIT NONE
  INTEGER,PARAMETER :: null = 999
  REAL(4),PARAMETER :: Undef = 999
contains

! subroutine get_inverse_matrix(G, G_inv, N)
!   implicit none
!   integer, intent(in) :: N
!   double precision, intent(in) :: G(N, N)
!   double precision, intent(out) :: G_inv(N, N)
!   double precision,allocatable :: work(:)
!   integer :: ipiv(N), info, lwork
!   double precision :: G_temp(N, N)
!   G_temp = G
!   call dgetrf(N, N, G_temp, N, ipiv, info)
!   if (info /= 0) then
!     print *, "Error: dgetrf failed with info =", info
!     stop
!   end if
!   lwork = N
!   call dgetri(N, G_temp, N, ipiv, work, lwork, info)
!   if (info /= 0) then
!     print *, "Error: dgetri workspace query failed with info =", info
!     stop
!   end if
!   lwork = int(work(1))
!   allocate(work(lwork))
!   call dgetri(N, G_temp, N, ipiv, work, lwork, info)
!   if (info /= 0) then
!     print *, "Error: dgetri failed with info =", info
!     stop
!   end if
!   G_inv = G_temp
!   deallocate(work)
! end subroutine get_inverse_matrix


subroutine get_inverse_matrix(G,Ginv,info,n)
implicit none
integer,intent(in) :: n
real(8),intent(in) :: G(n,n)
real(8),intent(out) :: Ginv(n,n)
integer,intent(out) :: info
integer :: lda
integer :: i,j
integer :: ipiv(n)
integer :: lwork
real(8) :: G2(n,n)
real(8),allocatable :: work(:)
  lda=n
  lwork=n
  Ginv = 0.0
  Ginv=G
  allocate(work(lwork))
  call dgetrf(n,n,Ginv,lda,ipiv,info)
  call dgetri(n,Ginv,lda,ipiv,work,lwork,info)
  !Ginv = G2
endsubroutine get_inverse_matrix



  SUBROUTINE matrix_product(A, B, C, M, N, K)
    implicit none
    REAL(8), INTENT(IN) :: A(M,K) 
    REAL(8), INTENT(IN) :: B(K,N) 
    REAL(8), INTENT(OUT) :: C(M,N) 
    INTEGER, INTENT(IN) :: M, N, K 

    INTEGER :: i, j, l

    C = 0.0D0

    DO i = 1, M
      DO j = 1, N
        DO l = 1, K
          C(i, j) = C(i, j) + A(i, l) * B(l, j)
        END DO
      END DO
    END DO
  END SUBROUTINE matrix_product

  SUBROUTINE matrix_vector_product(A, x, y, M, N)
  implicit none
  integer, intent(in) :: M, N
  REAL(8), intent(in) :: A(M, N), x(N)
  REAL(8), intent(out) :: y(M)
  INTEGER :: i, j, l
    y = 0.0D0
    ! 行列積の計算
    DO i = 1, M
      DO j = 1, N
          y(i) = y(i) + A(i, j) * x(j)
      END DO
    END DO
  END SUBROUTINE matrix_vector_product

! subroutine matrix_product( A, B, C,M, N, K)
!   implicit none
!   integer, intent(in) :: M, N, K
!   double precision, intent(in) :: A(M, K), B(K, N)
!   double precision, intent(out) :: C(M, N)
!   double precision :: alpha, beta
!   integer :: lda, ldb, ldc
!   alpha = 1.0d0
!   beta = 0.0d0
!   lda = M
!   ldb = K
!   ldc = M
!   call dgemm('N', 'N', M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)
! end subroutine matrix_product

! subroutine matrix_vector_product(A, x, y, M, N)
!   implicit none
!   integer, intent(in) :: M, N
!   double precision, intent(in) :: A(M, N), x(N)
!   double precision, intent(out) :: y(M)
!   double precision :: alpha, beta
!   integer :: lda
!   double precision :: work(M)
!   alpha = 1.0d0  ! Scaling factor for A * x
!   beta = 0.0d0   ! Scaling factor for y (initially zero)
!   lda = M        ! Leading dimension of A
!   call dgemv('N', M, N, alpha, A, lda, x, 1, beta, y, 1)
! end subroutine matrix_vector_product

! compute eigenvalue & eigenvector by Jabobi method
subroutine get_eigen_value_by_Jacobi(A,X,n)
! state variables
implicit none
integer :: i,k,l
integer, intent(in) :: n
integer :: p,q,fg
real(8) :: bg,ss,tt,vv,sai,sg,csai
real(8) :: bpi,bqi,bpp,bqq,xpq,xqp
real(8) :: eps
real(8) :: A(n,n)
real(8) :: X(n,n)
! n is the degree of matrix (You should check.)
! allocate(A(n,n))
! allocate(X(n,n))

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

end subroutine get_eigen_value_by_Jacobi

SUBROUTINE Inner_Product(inval,isize2,jsize2,ksize2,inlat,inlon,isize,jsize,ymin,xmin,yint,xint,output)
   IMPLICIT NONE
   INTEGER :: i,j,ii,jj,k
   iNTEGER :: ipos,jpos,ipos2,jpos2
   REAL(4) :: dx,dy
   INTEGER(4),INTENT(IN) :: jsize,isize,isize2,jsize2,ksize2
   REAL(4),INTENT(IN) :: xmin,ymin,xint,yint
   REAL(4),INTENT(IN) :: inlat(1:isize,1:jsize),inlon(1:isize,1:jsize)
   REAL(4),INTENT(IN) :: inval(1:ksize2,1:isize2,1:jsize2)
   REAL(4) :: rlon(1:jsize2),rlat(1:isize2)
   REAL(4),INTENT(OUT):: output(1:ksize2,1:isize,1:jsize)

   rlon(1:jsize2) = (/ (xmin+xint*(i-1),i=1,jsize2) /)
   where(rlon > 360.)
      rlon=rlon-360.
   endwhere
   rlat(1:isize2) = (/ (ymin+yint*(i-1),i=1,isize2) /)
   print *,isize,jsize,isize2,jsize2,xmin,ymin,xint,yint
   
   do j = 1, jsize
      do i = 1, isize
         jpos=-999
         do jj = 2, jsize2
            if((inlon(i,j) > rlon(jj-1)) .and. &
                  (inlon(i,j) <= rlon(jj))) then
               jpos=jj
               jpos2=jj-1
            endif
         enddo


      enddo
    enddo
    
END SUBROUTINE Inner_Product 

function bilin(dx,dy,Undef,var1,var2,var3,var4)
  implicit none
  REAL(4) :: bilin
  REAL(4),INTENT(in) :: dx, dy, Undef
  REAL(4),INTENT(in) :: var1,var2,var3,var4
!   REAL(4),INTENT(in) :: la1,la2,la3,la4

  if((var1 /= Undef) .and. (var2 /= Undef) .and. &
       & (var3 /= Undef) .and. (var4 /= Undef)) then
    bilin = (1.-dx)*(1.-dy)*var1 &
         & +    dx *(1.-dy)*var2 &
         & +    dx *    dy *var3 &
         & +(1.-dx)*    dy *var4
  else
    bilin = Undef
  endif
end function bilin


subroutine deming_linear_from_buf(x, y, sx, sy, n, b0, b1)
  ! 線形 EIV（Deming）の閉形式解。delta = mean(sy^2)/mean(sx^2)
  implicit none
  integer(4), intent(in) :: n
  real(4),    intent(in) :: x(n), y(n), sx(n), sy(n)
  real(8),    intent(out):: b0, b1
  integer(4) :: i
  real(8) :: Ex, Ey, Sxx, Syy, Sxy, vx, vy, delta, D, sxy_eps

  if (n < 2) then
    b0 = 0.0d0; b1 = 0.0d0; return
  end if

  Ex = 0.0d0; Ey = 0.0d0; vx = 0.0d0; vy = 0.0d0
  do i = 1, n
    Ex = Ex + dble(x(i)); Ey = Ey + dble(y(i))
    vx = vx + dble(sx(i))**2; vy = vy + dble(sy(i))**2
  end do
  Ex = Ex/dble(n); Ey = Ey/dble(n)
  vx = vx/dble(n); vy = vy/dble(n)
  delta = vy / (vx + 1.0d-12)

  Sxx = 0.0d0; Syy = 0.0d0; Sxy = 0.0d0
  do i = 1, n
    Sxx = Sxx + (dble(x(i))-Ex)**2
    Syy = Syy + (dble(y(i))-Ey)**2
    Sxy = Sxy + (dble(x(i))-Ex)*(dble(y(i))-Ey)
  end do
  Sxx = Sxx/dble(n); Syy = Syy/dble(n); Sxy = Sxy/dble(n)

  sxy_eps = max(1.0d-12, abs(Sxy))
  D  = Syy - delta*Sxx
  b1 = ( D + sqrt( max(0.0d0, D*D + 4.0d0*delta*Sxy*Sxy) ) ) / ( 2.0d0 * sxy_eps ) * sign(1.0d0, Sxy)
  b0 = Ey - b1*Ex
end subroutine deming_linear_from_buf

subroutine solve3x3_sym(H, b, x)
  ! 対称 3×3 の連立をガウス消去で解く（簡易）
  implicit none
  real(8), intent(inout) :: H(3,3)
  real(8), intent(in)    :: b(3)
  real(8), intent(out)   :: x(3)
  real(8) :: A(3,3), rhs(3), piv, m
  integer :: i, j, k

  A = H; rhs = b

  do k = 1, 3
    piv = A(k,k); if (abs(piv) < 1.0d-18) piv = sign(1.0d-18, piv)
    do i = k+1, 3
      m = A(i,k)/piv
      do j = k, 3
        A(i,j) = A(i,j) - m*A(k,j)
      end do
      rhs(i) = rhs(i) - m*rhs(k)
    end do
  end do

  do i = 3, 1, -1
    piv = A(i,i); if (abs(piv) < 1.0d-18) piv = sign(1.0d-18, piv)
    x(i) = (rhs(i) - sum( A(i,i+1:3)*x(i+1:3) )) / piv
  end do
end subroutine solve3x3_sym



end module
