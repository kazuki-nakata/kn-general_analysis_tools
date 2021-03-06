Program BBS
    IMPLICIT NONE
  INTEGER :: i,j,ii,jj,k,kk,iii,jjj,i2,j2,i3,j3,i4,j4,scan_i,scan_j
  INTEGER(4) :: ii1,ii2,jj1,jj2,kx,ky,window
  INTEGER(4) :: inum
  INTEGER(4) :: mi,mj,isize2,jsize2,rad,nt,step      
  INTEGER(4),PARAMETER :: isize=100
  INTEGER(4),PARAMETER :: jsize=100
  INTEGER(4),PARAMETER :: band=4
  INTEGER(4),PARAMETER :: imax_tl=5
  INTEGER(4),PARAMETER :: jmax_tl=5
  INTEGER(4),PARAMETER :: n_templ=1
  REAL(4) :: templ(imax_tl,jmax_tl,band,n_templ)
  INTEGER(4) :: minx2(imax_tl,jmax_tl), miny2(imax_tl,jmax_tl)  
  REAL(4) :: data_i(isize,jsize,band)
  REAL(4) :: oscale(isize,jsize)
  REAL(4) :: orient(isize,jsize)
  REAL(4) :: tmp_orient(n_templ)
  INTEGER(4) :: mask(isize,jsize)
  INTEGER(4),ALLOCATABLE :: minx(:,:), miny(:,:)  
  REAL(4),ALLOCATABLE :: x1(:,:),x2(:,:),x3(:,:,:),x4(:,:,:)
  REAL(4),ALLOCATABLE :: y1(:,:),y2(:,:),y3(:,:,:),y4(:,:,:)
  real(4),ALLOCATABLE :: ndvi1(:,:,:),si1(:,:,:),conv(:,:,:)
  REAL(4),PARAMETER :: disdep= 8.0  !dependency factor of pixels distance
  INTEGER(4),PARAMETER :: psize=1
  REAL(4), PARAMETER :: cal_n=15
  REAL(4) :: ai1,ai2,aj1,aj2,dis_rgb,pi,dora,theta
  REAL(4) :: dis_xy(imax_tl,jmax_tl)
  REAL(4) :: distance,maxdis,adis,bbp_max,asphalt
  real(4) :: dis(imax_tl,jmax_tl)
  real(4) :: sum,arad
  real(4) :: ndvi2,si2
  integer(4) :: tespos
  tespos=1
  templ(:,:,:,:)=5
  tmp_orient(:)=0
  oscale(:,:)=2
  data_i(:,:,:)=5
  mask(:,:)=0
  orient(:,:)=45
  pi=acos(-1.)
  dora=pi/180.

  kx=(imax_tl-1)/2
  ky=(jmax_tl-1)/2
  theta=0
  
  allocate(x3(imax_tl,jmax_tl,n_templ),x4(imax_tl,jmax_tl,n_templ))
  allocate(y3(imax_tl,jmax_tl,n_templ),y4(imax_tl,jmax_tl,n_templ))
  allocate(ndvi1(imax_tl,jmax_tl,n_templ),si1(imax_tl,jmax_tl,n_templ))


  do nt=1,n_templ
  maxdis=sqrt(real(kx)**2+real(ky)**2)  
        do jj = 1, jmax_tl
           do ii = 1, imax_tl
              i2=ii-kx-1
              j2=jj-ky-1
              ai2=real(ii)-real(kx)-1
              aj2=real(jj)-real(ky)-1
           x3(ii,jj,nt)=ai2*cos(theta)+aj2*sin(theta)
           y3(ii,jj,nt)=aj2*cos(theta)-ai2*sin(theta)
           x4(ii,jj,nt)=(ai2*cos(theta)+aj2*sin(theta))/real(imax_tl-1)
           y4(ii,jj,nt)=(aj2*cos(theta)-ai2*sin(theta))/real(jmax_tl-1)
           ndvi1(ii,jj,nt)=(templ(ii,jj,4,nt)-templ(ii,jj,3,nt)) &
                        /(templ(ii,jj,4,nt)+templ(ii,jj,3,nt))
           si1(ii,jj,nt)=atan(templ(ii,jj,1,nt) &
                /max(templ(ii,jj,2,nt),templ(ii,jj,3,nt)))+0.25
        enddo
        print *,x4(:,jj,1)
!        print *, '-------------------'
!        print *,y4(:,jj,1)
       enddo
   enddo



   
    do j =50, 50
      do i =50, 50

         window=nint(oscale(i,j)*2.0*2.0)   !1.5 times radius
         step=int(real(window)/cal_n)


         if(step.le.0)  step=1
         window=((int((window-1)/step+1)-1)*step+1)

         if(step.eq.1) then
            mi=mod(window,2)
            if(mi.eq.0) window=window+1
         endif

         arad=(real(window-1))/2


         allocate(x1(window,window),y1(window,window))
         allocate(x2(window,window),y2(window,window))
         allocate(minx(window,window),miny(window,window))
         allocate(conv(window,window,4))
         print *,window,step,arad
         do scan_j=2,2
            do scan_i=2,2                  
               i4=i+int(arad/4.)*(scan_i-2)
               j4=j+int(arad/4.)*(scan_j-2)
               print *,i4,j4
               if((i4.lt.1).or.(i4.gt.isize).or.(j4.lt.1).or.(j4.gt.jsize)) cycle

               bbp_max=0
               do nt =1, n_templ
                  do kk =1,1        
                     if(kk.eq.1) theta=(tmp_orient(nt)-orient(i,j))*dora
                     if(kk.eq.2) theta=((tmp_orient(nt)-orient(i,j))-180)*dora
                     !print *,'---------------------------------------------'
                     print *,theta
                     do jj = 1, window
                        do ii = 1, window
                           ai2=real(ii)-arad-1
                           aj2=real(jj)-arad-1
                           x1(ii,jj)=ai2*cos(theta)+aj2*sin(theta)
                           y1(ii,jj)=aj2*cos(theta)-ai2*sin(theta)
                           x2(ii,jj)=(ai2*cos(theta)+aj2*sin(theta))/(real(window)-1)
                           y2(ii,jj)=(aj2*cos(theta)-ai2*sin(theta))/(real(window)-1)
                           !              print *,x1(ii,jj),x2(ii,jj),ii,jj
                        enddo
                        print *,x2(:,jj)
                        
                     enddo
!!$                    ! print *,'---------------------------------------------'
!!$                     do jj = 1, window
!!$                        print *,y2(:,jj)
!!$                     enddo
                     
                     conv(:,:,:)=0
                     do jj=1,window,step
                        do ii=1,window,step
 !                          print *,ii1,jj1
                           ii1=int(x1(ii,jj))+i4;jj1=int(y1(ii,jj))+j4
                           sum=0
                           do jjj =1,step
                              do iii =1,step
                                 ii2=(ii1-1)+iii-int(real(step)/2)
                                 jj2=(jj1-1)+jjj-int(real(step)/2)
                                 if((ii2.gt.isize).or.(ii2.lt.1)) cycle
                                 if((jj2.gt.jsize).or.(jj2.lt.1)) cycle
                                 sum=1+sum
                                 do k =1,4
                                    conv(ii,jj,k)=conv(ii,jj,k)+data_i(ii2,jj2,k)
                                 enddo
                              enddo
                           enddo
                           do k =1,4
                              conv(ii,jj,k)=conv(ii,jj,k)/sum
                              !                 print *,conv(ii,jj,k)
                           enddo
                        enddo
                     enddo

                     !----------------------estimate BBPs-------------------------------
                     inum=0
                     minx(:,:)=999
                     miny(:,:)=999
                     dis(:,:)=10000.
                     sum=0
                     !        print *,i,j,step,window

                     do jj=tespos,tespos,step
                        do ii=tespos,tespos,step
                           distance=1000000.
                           ii1=int(x1(ii,jj))+i4;jj1=int(y1(ii,jj))+j4
                           print *,ii,jj,ii1,jj1,1
                           ai1=x2(ii,jj);aj1=y2(ii,jj)
                           print *,ai1,aj1
                           if((ii1.lt.1).or.(ii1.gt.isize).or.(jj1.lt.1).or.(jj1.gt.jsize)) cycle
                           if((ai1.lt.-0.51).or.(ai1.gt.0.51).or.(aj1.lt.-0.51).or.(aj1.gt.0.51)) cycle
                           sum=sum+1              
                           do jjj=1,jmax_tl,psize
                              do iii=1,imax_tl,psize
                                 ii2=int(x3(iii,jjj,nt));jj2=int(y3(iii,jjj,nt))
                                 ai2=x4(iii,jjj,nt);aj2=y4(iii,jjj,nt)
                                 !                                 print *,ai1,ai2,aj1,aj2
                                 
                                 dis_xy(iii,jjj)=sqrt(real(ai1-ai2)**2+real(aj1-aj2)**2)
                                 dis_rgb=0
                                 do k=1,band
                                    dis_rgb=(conv(ii,jj,k)-templ(iii,jjj,k,nt))**2 + dis_rgb
                                    !                       dis_rgb=(conv(ii,jj,k)/255.-templ(iii,jjj,k,nt)/255.)**2 + dis_rgb                      
                                 enddo

!                                 adis=dis_rgb+dis_xy*disdep

                                 if(adis.le.dis(iii,jjj)) then
                                    minx2(iii,jjj)=ii
                                    miny2(iii,jjj)=jj
                                    dis(iii,jjj)=adis
                                 endif
                                 if(adis.le.distance) then
                                    minx(ii,jj) = iii
                                    miny(ii,jj) = jjj
                                    distance=adis
                                 endif
                              enddo
                              print *,dis_xy(:,jjj)
                           enddo

                        enddo
                     enddo




                  enddo
               enddo



            enddo
         enddo

         deallocate(x1,y1,x2,y2,minx,miny,conv)

      enddo
   enddo  

 end Program BBS
 
SUBROUTINE Gaussian_kernel(output,wsize,sigma,sum,sf)
  IMPLICIT NONE
  INTEGER :: i,j,k
  REAL(4) :: pi,alpha,x,y,sigma2,sum,sf
  REAL(4),INTENT(IN) :: sigma
  INTEGER(4),INTENT(IN) :: wsize
  REAL(4),INTENT(OUT):: output(1:wsize,1:wsize)
  pi=acos(-1.)
  sigma2=2.*sigma**2
  alpha=1/(pi*sigma2)
  sum=0
  do j =1,wsize
     do i =1,wsize
        x=real(i)-(wsize+1)/2.
        y=real(j)-(wsize+1)/2.
        output(i,j)=alpha*exp(-((x*sf)**2+(y*sf)**2)/(sigma2))
        sum=sum+output(i,j)
     enddo
  enddo
!  print *,alpha,sum
END SUBROUTINE Gaussian_kernel
