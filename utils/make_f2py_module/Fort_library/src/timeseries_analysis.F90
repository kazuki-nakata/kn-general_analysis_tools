MODULE timeseries_analysis
  USE Basic_processing
  use OBJECT_MODULE3
  IMPLICIT NONE
!  INTEGER, PARAMETER :: MAX_IO = 64

CONTAINS

  SUBROUTINE incbeta_openland_detection(data_i,data_o,prior_beta1,prior_beta2,band,jsize,isize,p1,p2,num_pre,num_post)
  IMPLICIT NONE
  INTEGER :: i,j,k,t
  INTEGER :: ii,jj,kk
  INTEGER(4),INTENT(IN) :: isize,jsize,band
  REAL(8), INTENT(IN) :: p1,p2 !p1 and p2 are prob thresholds
  INTEGER(4), INTENT(OUT)  :: data_o(3,jsize,isize)
  INTEGER(1), INTENT(IN) :: data_i(band,jsize,isize)
  REAL(4), INTENT(IN) :: prior_beta1(2,jsize,isize)
  REAL(4), INTENT(IN) :: prior_beta2(2,jsize,isize)  
  INTEGER(4) :: probmap(jsize,isize)
  INTEGER(4) :: num_pre,num_post
  INTEGER(1) :: ts(band),ts_index(band)
  INTEGER :: inum,irec,ui
  INTEGER :: sum
  REAL(8) :: prob1,prob2,prob_t
  REAL(8) :: param_a1,param_b1
  REAL(8) :: param_a2,param_b2  
  REAL(8) :: m,n
  REAL(8) :: size,prob_ave
  REAL(8) :: max_prob_t,argmax_prob1,argmax_prob2,argmax_t,openland,min_kld,kld,q
  INTEGER :: argmax_threshold
  REAL(8) :: sc
!  REAL(8) :: gamma_d,hygeo,incbeta
  
  integer(4), allocatable :: pairs1(:),pairs2(:)
  integer(4) :: ipos_8c(4),jpos_8c(4)  
  INTEGER(4) :: irec1,irec2
  INTEGER(4) :: dummy1(1000000),dummy2(1000000)
  integer(4) :: pmax_uf

  print *, prior_beta1(1:2,1,1),band,jsize,isize,p1,p2,num_pre,num_post
  data_o(:,:,:)=0
  
!  sc=0.00001
  
  do i = 1, isize
     print *,i
     do j = 1, jsize

        
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
        min_kld=10000.
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
           
           prob1=1-incbeta(m,n,p1)
           
           sum=0
           do k =t+1,inum
              if(ts(k).eq.1) sum=sum+1
           enddo
           
           m=dble(sum)+param_a2
           n=dble(inum)-dble(t)-dble(sum)+param_b2
           
           prob2=1-incbeta(m,n,p2)
           
           prob_t = prob1 * prob2
           q=(prob2+prob1)/2.           
!           kld = prob2*(log(prob2)-log(q))+(1-prob2)*(log(1-prob2)-log(1-q)) &
!                + prob1*(log(prob1)-log(q))+(1-prob1)*(log(1-prob1)-log(1-q))
!           if(kld.lt.min_kld) then   
           if(prob_t.ge.max_prob_t) then
              argmax_t= ts_index(t + 1)
              max_prob_t = prob_t
              argmax_prob1 = prob1
              argmax_prob2 = prob2
           endif

        enddo

        if((argmax_t.ge.num_pre).and.(argmax_t.le.num_post) &
             .and.(max_prob_t.gt.0.5)) then
           probmap(j,i)=int(max_prob_t*100)
        else
           probmap(j,i)=0
        endif

        data_o(1,j,i)=int(100*max_prob_t)
        data_o(2,j,i)=argmax_t

        
     enddo
  enddo


!---------------------Area average for probability-------------  
  call init_unionfind(jsize,isize,probmap,1,0,smax_uf)
  ipos_8c=(/-1,0,1,1/)
  jpos_8c=(/1,1,1,0/)
  allocate(pairs1(smax_uf*4),pairs2(smax_uf*4))
  pairs1(:)=-999;pairs2(:)=-999

  print *, 'make pixel pairs'
  irec=0
  do j=1,isize
     do i=1,jsize
        !--------------------8-connected system-------------------------
        if(ptrt(i,j).eq.0) cycle
        do kk=1,4
           if((i+ipos_8c(kk).le.0).or.(i+ipos_8c(kk).gt.isize)) cycle
           if((j+jpos_8c(kk).le.0).or.(j+jpos_8c(kk).gt.jsize)) cycle           
           if(ptrt(i+ipos_8c(kk),j+jpos_8c(kk)).ne.0) then
              irec=1+irec
              pairs1(irec)=ptrt(i,j)
              pairs2(irec)=ptrt(i+ipos_8c(kk),j+jpos_8c(kk))
           endif
        enddo

     enddo
  enddo

  pmax_uf=irec
  
  do i = 1,pmax_uf
     call union_uf(pairs1(i),pairs2(i),ui)
  enddo
  
!  data_o(3,:,:)=0
  do j=1,isize
     do i=1,jsize
        if(probmap(i,j).eq.0) cycle
        ii=find_uf(ptrt(i,j))
        size = get_size_uf(ii)
        prob_ave = get_property_uf(ii,1)
        if(size * prob_ave .ge. 80) data_o(3,i,j)=prob_ave
     enddo
  enddo


  call deallocate_uf
  deallocate(pairs1,pairs2)
  
END SUBROUTINE incbeta_openland_detection  

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

!!$REAL(8) FUNCTION gamma_d(a)
!!$  implicit none
!!$  REAL(8),intent(in) :: a
!!$  INTEGER :: k
!!$  INTEGER :: i,j
!!$  real(8) :: pi
!!$  real(8),PARAMETER:: e=2.71828182846
!!$  real(8) :: a2
!!$
!!$  pi=dble(acos(-1.))
!!$
!!$  if(a.ge.1) then
!!$  gamma_d = sqrt(2*pi/a)*(a/e)**a*(1.+1./12./a+1./288./a/a)
!!$else
!!$   a2=a+1.0
!!$   gamma_d = sqrt(2*pi/a2)*(a2/e)**a2*(1.+1./12./a2+1./288./a2/a2)/a
!!$   endif
!!$!  print *,gamma_d
!!$  
!!$end FUNCTION gamma_d

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
!  REAL(8) :: gamma_d,hygeo
  
  x1 = hygeo(dble(1.),m+n,m+1,q)
  x2 = (q**m)*((1-q)**n)/m
  x3=gamma_d(m+n)/gamma_d(m)/gamma_d(n)
  incbeta=x1*x2*x3

end FUNCTION incbeta

END MODULE timeseries_analysis
