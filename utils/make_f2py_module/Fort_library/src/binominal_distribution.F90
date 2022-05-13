program simplified_switching_model
  implicit none
  integer :: i,j
  integer :: num
!  integer :: data(num)
  integer :: threshold
  integer :: sum
  real(8) :: prob1,prob2,prob_t
  real(8) :: prob_a1,prob_a2
  real(8),parameter :: p1=0.8
  real(8),parameter :: p2=0.2
  real(8) :: p3
  real(8) :: q1,q2,param_a,param_b
  integer :: m,n
  real(8) :: x1,x2,x3
  real(8) :: max_prob
  integer :: argmax_threshold
  real(8) :: gamma_d,hygeo

  num=1
  sum=1
  q1=1.0
  q2=0.5
  param_a=1
  param_b=1
  m=int(sum+param_a-1)
  n=int(num-sum+param_b-1)

  print *,m,n

!  do i =1,100
  do j =0,10
     q2=0.1*real(j)
     prob1=0
  
  x1 = hygeo(dble(1.),dble(m+n),dble(m+1),q2)
  x2 = (q2**m)*((1-q2)**n)/m
  prob1=x1*x2
  x3=gamma_d(m+n)/gamma_d(m)/gamma_d(n)
  prob1=prob1*x3
  print *,q2,prob1,x1,x2,x3
  enddo
!  data(1:2)=0;data(4:6)=1!;data(15)=1;data(16)=1 
!  print *,data(:)
!  enddo
!!$  do j =1, n-1
!!$
!!$  threshold=j
!!$  !-------split data---------------
!!$  sum=0
!!$  do i =1,threshold
!!$     if(data(i).eq.1) sum=sum+1
!!$  enddo
!!$  call binominal_distribution(p1,threshold,sum,prob_a1,1)
!!$  call binominal_distribution(p2,threshold,sum,prob_a2,1)
!!$!  print *,p2,threshold,threshold-sum,sum
!!$
!!$  do i = 1,11
!!$     p3 = 0.1 * real(i - 1)
!!$     call binominal_distribution(p3,threshold,sum,prob_a1,1)
!!$     print *,p3,prob_a1
!!$  enddo
!!$  
!!$!  prob1=log(prob_a1)-log(prob_a1+prob_a2)
!!$  prob1=prob_a1/(prob_a1+prob_a2)
!!$!  print *, threshold,prob_a1,prob_a2,prob1,sum
!!$!  print *,sum
!!$!  print *, "second time series"
!!$  sum=0
!!$  do i =threshold+1,n
!!$     if(data(i).eq.1) sum=sum+1
!!$  enddo
!!$  !  print *,sum
!!$!  print *,i
!!$
!!$     
!!$  call binominal_distribution(p2,n-threshold,sum,prob_a2,1)
!!$  call binominal_distribution(p1,n-threshold,sum,prob_a1,1)
!!$!  print *,prob2
!!$
!!$  prob2=prob_a2/(prob_a1+prob_a2)
!!$! print *, threshold,prob_a1,prob_a2,prob2,sum
!!$  
!!$  prob_t=prob1 * prob2 !(prob2+prob1-1)/prob1
!!$!   print *, j,prob1,prob2,prob_t
!!$
!!$  if(prob_t.ge.max_prob) then
!!$     argmax_threshold=threshold
!!$     max_prob = prob_t
!!$  endif
!!$
!!$  enddo
!!$!----------------------------------
!!$  print *,max_prob,argmax_threshold

  
endprogram simplified_switching_model

subroutine binominal_distribution(p,n,n_thres1,prob,mode)
  implicit none
  real(8),intent(in) :: p
  integer,intent(in) :: n
  integer :: k
  integer :: i,j
  real(8) :: bin(0:n)
  real(8) :: n1,n2
  integer :: nmin,nmax
  integer :: mode
  real,parameter :: threshold=0.9
  integer,intent(in) :: n_thres1
  real(8),INTENT(OUT) :: prob
  
  do k = 0,n
     if(k.eq.0) then
        n1=1
        n2=1
     else
        n1=1
        do j =n-k+1,n
           n1=n1*real(j)
        enddo
        n2=1
        do j =1,k
           n2=n2*real(j)
        enddo
     endif
     
     bin(k)=n1/n2*(p**k)*(1-p)**(n-k)
  enddo

  if(mode.eq.1) then
     nmin = 0
     nmax = n_thres1
  else
     nmin = n_thres1
     nmax = n
  endif
  
  prob=0
  do k =nmin,nmax
     prob=bin(k)+prob
  enddo
  prob=bin(n_thres1)
endsubroutine binominal_distribution

REAL(8) FUNCTION gamma_d(a)
  implicit none
  integer,intent(in) :: a
  integer :: k
  integer :: i,j
  gamma_d = 1
  do i = 1, a-1
     gamma_d = dble(i) * gamma_d
  enddo
end FUNCTION gamma_d

REAL FUNCTION beta(a, b)
  implicit none
  integer,intent(in) :: a
  integer,intent(in) :: b
  real :: gamma_x,gamma_y,gamma_xy
  integer :: k
  integer :: i,j

  gamma_x = 1 
  do i = 1, a-1
     gamma_x = real(i) * gamma_x
  enddo
  
  gamma_y = 1 
  do i = 1, b-1
     gamma_y = real(i) * gamma_y
  enddo
  
  gamma_xy = 1 
  do i = 1, a+b-1
     gamma_xy = real(i) * gamma_xy
  enddo

  beta=gamma_x*gamma_y/gamma_xy
  
end FUNCTION beta

REAL(8) FUNCTION hygeo(a, b, c, x)
  implicit none
  real(8),intent(in) :: a,b,c,x
  integer :: k
  integer :: i,j
  real(8) :: p1,f,p0
  f = 0
  p0 = 1
  do i = 0, 1000
     p1=(a+dble(i))*(b+dble(i))/(c+dble(i))/(1+dble(i))*x*p0
     f = f + p0
     p0=p1
     if(p0.le.10E-3) then
        print *,i
        exit
     endif
  enddo
  hygeo=f
end FUNCTION hygeo

!!$ REAL(8) FUNCTION incbeta(m, n, q)
!!$  implicit none
!!$  real(8),intent(in) :: m,n,q
!!$  integer :: k
!!$  integer :: i,j
!!$  real(8) :: p1,f,p0
!!$  real(8) :: x1,x2,x3
!!$  
!!$  x1 = hygeo(dble(1.),dble(m+n),dble(m+1),q)
!!$  x2 = (q2**m)*((1-q2)**n)/m
!!$  x3=gamma_d(m+n)/gamma_d(m)/gamma_d(n)
!!$  incbeta=x1*x2*x3

!end FUNCTION incbeta
