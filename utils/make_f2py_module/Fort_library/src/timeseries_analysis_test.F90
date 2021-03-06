Program incbeta_openland_detection
  IMPLICIT NONE
  INTEGER :: i,j,k,t
  INTEGER(4),parameter :: isize=24
  INTEGER(4),PARAMETER :: jsize=24
  INTEGER(4),PARAMETER :: band=10
  REAL(8), PARAMETER :: p1 =0.8
  REAL(8), PARAMETER :: p2 =0.8
  !p1 and p2 are prob thresholds
  INTEGER(4) :: aaa
!  data aaa/0, 0, 0, 1, 1, 1, 1/
  REAL(4)  :: data_o(5,jsize,isize)
  INTEGER(1) :: data_i(band,jsize,isize)
  REAL(4) :: prior_beta1(2,jsize,isize)
  REAL(4) :: prior_beta2(2,jsize,isize)  
  INTEGER(1) :: ts(band),ts_index(band)
  INTEGER :: inum
  INTEGER :: sum
  REAL(8) :: prob1,prob2,prob_t
  REAL(8) :: param_a1,param_b1
  REAL(8) :: param_a2,param_b2
  REAL(8) :: m,n
  REAL(8) :: max_prob_t,argmax_prob1,argmax_prob2,argmax_t,openland,min_kl,kld,q
  INTEGER :: argmax_threshold
  REAL(8) :: gamma_d,hygeo,incbeta
  REAL(8) :: sc

!!$  do i = 1, band
!!$     data_i(i,:,:)=a(i)
!!$  enddo
  data_i(1:1,:,:)=0
  data_i(2:7,:,:)=1
  data_i(8:10,:,:)=1
!  data_i(3:5,:,:)=1
!  data_i(6:6,:,:)=0
!  data_i(7:7,:,:)=1
!  data_i(4,:,:)=1
!  data_i(9:12,:,:)=1
!  data_i(13:16,:,:)=1

!  data_i(2,:,:)=1
!  data_i(5,:,:)=1
!!$  data_i(1:1,:,:)=0
!!$  data_i(2:4,:,:)=1
!!$  
!!$  data_i(1:5,:,:)=0
!!$  data_i(6:10,:,:)=1
!!$
!!$  
!!$  !--32 data sampling case-----------
!!$  data_i(1:16,:,:)=0  
!!$  data_i(17:32,:,:)=1
!!$!--
!!$  data_i(1:16,:,:)=0
!!$  data_i(17:32,:,:)=1
!!$  do i =1,16,4
!!$     data_i(i,:,:)=1
!!$  enddo
!!$  do i =16,32,4
!!$     data_i(i,:,:)=0
!!$  enddo  
  !---seasonal change--------
!!$  data_i(1:4,:,:)=0
!!$  data_i(5:8,:,:)=1
!!$  data_i(9:12,:,:)=0
!!$  data_i(13:16,:,:)=1
!!$  data_i(17:20,:,:)=0
!!$  data_i(21:24,:,:)=1
!!$  data_i(25:28,:,:)=0
!!$  data_i(29:32,:,:)=1
  !------------------------
  !---seasonal change--------
!!$  data_i(1:8,:,:)=1
!!$  data_i(9:16,:,:)=0
!!$  data_i(17:24,:,:)=1
!!$  data_i(25:32,:,:)=0
!  data_i(28,:,:)=0
!  data_i(1:1,:,:)=1
!  data_i(2:2,:,:)=1
!  data_i(10:10,:,:)=1
  prior_beta1(1,:,:)=5.5
  prior_beta1(2,:,:)=1/2.
  prior_beta2(1,:,:)=1/2.
  prior_beta2(2,:,:)=1/2.  
!  data_i(:,:,:)=1
!  data_o(:,:,:)=-1.
  open(30,file="analysis.txt",status="replace",form='formatted', &
       access='sequential')
  
  sc=0.0000
  do i = 1, 1
     do j = 1, 1
        
        inum=0
        do k =1, band
!           print *,data_i(k,j,i)
           if((data_i(k,j,i).eq.1).or.(data_i(k,j,i).eq.0)) then
              inum = 1 + inum              
              ts(inum) = data_i(k,j,i)
              ts_index(inum) = k
           endif
        enddo
        
        if(inum.le.1) cycle
        min_kl=10000.
        argmax_t= 0
        max_prob_t = 0
        argmax_prob1 = 0
        argmax_prob2 = 0
        param_a1 = dble(prior_beta1(1,j,i))
        param_b1 = dble(prior_beta1(2,j,i))
        param_a2 = dble(prior_beta2(1,j,i))
        param_b2 = dble(prior_beta2(2,j,i))
        
        do t = 0, inum-1

           sum=0
          
           do k =1, t
              if(ts(k).eq.0) sum=sum+1
           enddo

           
           m=dble(sum+param_a1)
           n=dble(t)-dble(sum)+param_b1
!           print *,m,n,'hi'
           prob1=1-incbeta(m,n,p1)
!           print *,m,n
!           print *,sum,t,j,prob1,m,n        
!           print *,sum,t,m,n,param_a,param_b
           
           sum=0
           do k =t+1,inum
              if(ts(k).eq.1) sum=sum+1
           enddo

           m=dble(sum)+param_a2
           n=dble(inum)-dble(t)-dble(sum)+param_b2
!           print *,m,n
           prob2=1-incbeta(m,n,p2)
           
           prob_t = prob1*prob2!prob2/abs(-log(prob2)+log(prob1))/100
           kld = prob2*(log(prob2)-log(prob1))+(1-prob2)*(log(1-prob2)-log(1-prob1))
           q=(prob2+prob1)/2.
           kld = prob2*(log(prob2) -log(q))+(1-prob2)*(log(1-prob2)-log(1-q)) &
                + prob1*(log(prob1)-log(q))+(1-prob1)*(log(1-prob1)-log(1-q))
 !          kld =(1-prob2)*(log(1-prob2)) &
 !               +(1-prob1)*(log(1-prob1))
!           kld = kld/3.
!           kld=-(log(prob1)+log(prob2))/10.
!           kld = prob2*log(prob2)+(1-prob2)*log(1-prob2) &
!                + prob1*log(prob1)+(1-prob1)*log(1-prob1)
!           kld = kld/4.
!           kld = prob2*log(prob2)/(prob1*log(prob1))
!           kld=log(kld)
!           prob_t = abs((-log(1-prob2)+log((prob1+1-prob2)/2.))*(1-prob2)/2. + prob_t)
           !kld=q
           print *,sum,t,j,data_i(t,1,1),prob_t,prob1,prob2,kld
           write(30,*) t,data_i(t,1,1),t+0.5,prob1,prob2,prob_t,kld
           if(prob_t.ge.max_prob_t) then
!           if(kld.lt.min_kl) then
              argmax_t= ts_index(t+1)
              max_prob_t = prob_t
              argmax_prob1 = prob1
              argmax_prob2 = prob2
              min_kl=kld
           endif

        enddo
        
        if(max_prob_t .ge.0.5) openland=1

        data_o(1,j,i) = openland
        data_o(2,j,i) = argmax_t
        data_o(3,j,i) = max_prob_t
        data_o(4,j,i) = argmax_prob1
        data_o(5,j,i) = argmax_prob2

     enddo
  enddo

  write(30,*) band,data_i(band,1,1)

  print *,data_o(:,1,1)
  
END PROGRAM incbeta_openland_detection  

subroutine binominal_distribution(p,n,n_thres1,prob,mode)
  implicit none
  REAL(8),intent(in) :: p
  INTEGER,intent(in) :: n
  INTEGER :: k
  INTEGER :: i,j
  REAL(8) :: bin(0:n)
  REAL(8) :: n1,n2
  INTEGER :: nmin,nmax
  INTEGER :: mode
  REAL,parameter :: threshold=0.9
  INTEGER,intent(in) :: n_thres1
  REAL(8),INTENT(OUT) :: prob
  
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
  REAL(8),intent(in) :: a
  INTEGER :: k
  INTEGER :: i,j
  real(8) :: pi
  real(8),PARAMETER:: e=2.718281828459045
  real(8) :: a2,b

  pi=dble(acos(-1.))

  if(a.ge.1) then
     gamma_d = sqrt(2*pi/a)*(a/e)**a
     b=(1.+1./12./a+1./288./a/a-139./51840./a/a/a-571./2488320./a/a/a/a)
     gamma_d =gamma_d*b
else
   a2=a+1.0
   gamma_d = sqrt(2*pi/a2)*(a2/e)**a2*(1.+1./12./a2+1./288./a2/a2-139./51840./a2/a2/a2-571./2488320./a/a/a/a)/a
   endif
!  print *,gamma_d
  
end FUNCTION gamma_d

REAL FUNCTION beta(a, b)
  implicit none
  INTEGER,intent(in) :: a
  INTEGER,intent(in) :: b
  REAL :: gamma_x,gamma_y,gamma_xy
  INTEGER :: k
  INTEGER :: i,j

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
  REAL(8),intent(in) :: a,b,c,x
  INTEGER :: k
  INTEGER :: i,j
  REAL(8) :: p1,f,p0
  f = 0
  p0 = 1
  do i = 0, 1000
     p1=(a+dble(i))*(b+dble(i))/(c+dble(i))/(1+dble(i))*x*p0
     f = f + p0
     p0=p1
     if(p0.le.10E-3) exit
  enddo
  hygeo=f
end FUNCTION hygeo

 REAL(8) FUNCTION incbeta(m, n, q)
  implicit none
  REAL(8),intent(IN) :: q
  REAL(8),INTENT(IN) :: m,n
  INTEGER :: k
  INTEGER :: i,j
  REAL(8) :: p1,f,p0
  REAL(8) :: x1,x2,x3
  REAL(8) :: gamma_d,hygeo
  
  x1 = hygeo(dble(1.),m+n,m+1,q)
  x2 = (q**(m))*((1-q)**(n))/m
  x3=gamma_d(m+n)/gamma_d(m)/gamma_d(n)
  incbeta=x1*x2*x3
!  print *,x1,x2,x3,m,n

end FUNCTION incbeta

