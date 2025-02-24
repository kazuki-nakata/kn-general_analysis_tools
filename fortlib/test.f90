program matinv


  implicit none

  real(8) :: M(4,4)
  real(8) :: Minv(4,4)
  integer :: lda =4
  integer :: info,i,j
  integer :: ipiv(4)
  integer :: lwork=4
  real(8),allocatable::work(:)
  integer :: n =4

  M(1,1) = 1.2
  M(1,2) = 2.3
  M(1,3) = 2.4
  M(1,4) = 3.2
  M(2,1) = 0.7
  M(2,2) = 0.23
  M(2,3) = 3.2
  M(2,4) = 1.2
  M(3,1) = 2.4
  M(3,2) = 0.1
  M(3,3) = 3.5
  M(3,4) = 0.1
  M(4,1) = 3.2
  M(4,2) = 0.3
  M(4,3) = -4.1
  M(4,4) = 0.0


  Minv = 0.0
  allocate (work(lwork))
  call dgetrf(n,n,M,lda,ipiv,info)
  call dgetri(n,M,lda,ipiv,work,lwork,info)

  Minv = M
  do i =1,4
     write(*,*) Minv(i,1:4)
  enddo
  end program matinv
