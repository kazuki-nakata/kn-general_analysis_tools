MODULE Template_Matching_lr
  USE Basic_Processing
  IMPLICIT NONE
  INTEGER, PARAMETER :: MAX_IO = 64

  CONTAINS

  SUBROUTINE landslide_detection(data_i,data_o,band,jsize,isize,threshold,tmp_i,num_t,band_t,jsize_t,isize_t)
  INTEGER :: i,j,k,l,ii,jj,kk,i2,j2,k2,i3,j3,k3      
  REAL, PARAMETER :: threshold_b= 0.3
  INTEGER, PARAMETER :: band1=3
  INTEGER, PARAMETER :: band2=4
  INTEGER(4),INTENT(IN) :: isize_t,jsize_t,band_t,num_t  
  INTEGER(4),INTENT(IN) :: isize,jsize,band
  INTEGER(4),INTENT(OUT) :: data_o(jsize,isize)
  REAL(4),INTENT(IN) :: data_i(band,jsize,isize)
  REAL(4), INTENT(IN) :: tmp_i(num_t,band_t,jsize_t,isize_t)  
  REAL(4),INTENT(IN) :: threshold
  INTEGER(4),ALLOCATABLE :: morp1(:,:),morp2(:,:),bin(:,:)
  REAL(4),ALLOCATABLE :: dum(:,:),ndvi(:,:),data(:,:,:)
  INTEGER(4),ALLOCATABLE :: mask1(:,:),mask2(:,:)
   REAL(4), ALLOCATABLE :: tmp0(:,:,:,:)
   REAL(4), ALLOCATABLE :: tmp_orient(:)
   INTEGER(4), PARAMETER :: imax_tl=17
   INTEGER(4), PARAMETER :: jmax_tl=17
   INTEGER(4) :: idirect,jdirect
   REAL(4), ALLOCATABLE :: tmp(:,:,:,:)  
!-----------DoG-------------------------------
  REAL(4),ALLOCATABLE :: octave(:,:,:),diff(:,:,:),oscale(:,:),orient(:,:)
  REAL(4),ALLOCATABLE :: dummy(:,:),confidence(:,:),bbp(:,:)
  INTEGER(4),ALLOCATABLE :: dummy2(:,:),forshape(:,:)
  INTEGER, PARAMETER :: octaves = 3
  INTEGER, PARAMETER :: interval = 2
  INTEGER, PARAMETER :: numpero = 6
  REAL, PARAMETER :: sigma=0.5
  REAL(4) :: sum,ksigma,k_int,diffmax,sf,diffmin
  INTEGER :: mi,mj,wsize,dum_scale,disize,djsize,n_oc,interval2
  INTEGER(4) :: isize2,jsize2,isize3,jsize3,conf_index,window
  character(512),INTENT(IN) :: dir

   !---------------INITIALIZATION-----------------------------
  allocate(ndvi(isize,jsize),bin(isize,jsize))
!!!  allocate(morp1(isize,jsize),morp2(isize,jsize))
  allocate(mask1(isize,jsize),mask2(isize,jsize),data(isize,jsize,band))
  CALL reshape_backforward_3D(data_i(1:band,1:jsize,1:isize),data(1:isize,1:jsize,1:band),band,jsize,isize)
  
  !--------------MAKE MASK DATA FOR SMALL SCALE REGION----------
  !---opening for removing small scale objects in high ndvi regions
  !---closing for removing small scale objects in low ndvi regions
  !---candidates of calculation estimation  
  
  CALL Index(data(1:isize,1:jsize,1:band),ndvi(1:isize,1:jsize),isize,jsize,band,band1,band2,0)
!  CALL Interpolation_type1(ndvi(1:isize,1:jsize),isize,jsize,999.)
  CALL Binarization(ndvi(1:isize,1:jsize),bin(1:isize,1:jsize),isize,jsize,threshold_b,999)
 
  CALL Make_Mask_FORLANDSLIDE(data(:,:,4:5),mask1(:,:),isize,jsize)
  CALL Mask_MERGE(bin(:,:),mask1(:,:),mask2(:,:),isize,jsize)
  !  morp1(:,:)=mask1(:,:)


  deallocate(bin)
  mask1(:,:)=0


!!$  CALL Index(data(:,:,:),ndvi(:,:),isize,jsize,band,band1,band2,0)
!!$!  CALL Interpolation_type1(ndvi(:,:),mask1(:,:),isize,jsize,999.)
!!$  CALL Binarization(ndvi(:,:),bin(:,:),isize,jsize,threshold_b,999)
!!$  CALL Morphology_Binary_Erosion(bin(:,:),morp1(:,:),isize,jsize,3,999)
!!$  CALL Morphology_Binary_Dilation(morp1(:,:),morp2(:,:),isize,jsize,3,999)
!!$  DO i =1,3
!!$     CALL Morphology_Binary_Dilation(morp2(:,:),morp1(:,:),isize,jsize,3,999)
!!$     morp2(:,:)=morp1(:,:)
!!$  ENDDO
!!$  DO i=1,3
!!$  CALL Morphology_Binary_Erosion(morp2(:,:),morp1(:,:),isize,jsize,3,999)
!!$  morp2(:,:)=morp1(:,:)
!!$  ENDDO
!!$  CALL Raster_Difference_Int(morp1(:,:),bin(:,:),mask2(:,:),isize,jsize,999)
!!$  deallocate(morp1,morp2,bin)
!!$
!-------------CREATE SCALE-SPACE MAPS BY DOG OPERATOR----------  

  wsize=6*interval+1
  k_int=(real(interval)/sigma)**(1/real(numpero))
  allocate(diff(isize,jsize,numpero*octaves))
  diff(:,:,:)=0.

  do n_oc=1, octaves ! octave number
     sf=2.**(real(n_oc)-1)
     interval2=interval**(n_oc-1)
     mi=mod(isize,interval2)
     mj=mod(jsize,interval2)
     isize2=(isize-mi)/interval2
     jsize2=(jsize-mj)/interval2
     allocate(octave(1:isize2,1:jsize2,0:numpero))
     if(n_oc.ne.1) allocate(mask1(1:isize2,1:jsize2))
     if(n_oc.eq.1) octave(:,:,0)=ndvi(:,:)
     if(n_oc.ne.1) CALL down_scaling_real(dummy(:,:),octave(:,:,0),disize,disize,isize2,jsize2,2)
     if(n_oc.ne.1) CALL down_scaling_int(dummy2(:,:),mask1(:,:),disize,disize,isize2,jsize2,2)
     do i =1, numpero
        if(n_oc.eq.1) ksigma = sigma*(k_int**real(i))
        if(n_oc.ne.1) ksigma= sqrt((real(interval2)*(k_int**real(i)))**2-real(interval2)**2)
        CALL gaussian_filter(octave(:,:,0),octave(:,:,i),mask1(:,:),isize2,jsize2,wsize,ksigma,sf) 
     enddo

     do k =1,numpero
        kk=k+(n_oc-1)*numpero
        do j =1,jsize
           do i =1,isize
              i2=(i-1)/interval2+1
              j2=(j-1)/interval2+1
              if((i2.gt.isize2).or.(j2.gt.jsize2)) cycle
              diff(i,j,kk)=octave(i2,j2,k)-octave(i2,j2,k-1)
           enddo
        enddo
     enddo

  !---------store for the next octave----------------------
     if(n_oc.ne.1) deallocate(dummy)
     if(n_oc.ne.1) deallocate(dummy2)  
     disize=isize2;djsize=jsize2
     allocate(dummy(disize,djsize))
     allocate(dummy2(disize,djsize))
     dummy(:,:)=octave(:,:,numpero)
     dummy2(:,:)=mask1(:,:)

     deallocate(octave,mask1)

  enddo
  
  !------scale selection by searching a maximum value----------
  allocate(confidence(isize,jsize))
  confidence(:,:)=0.
  allocate(oscale(isize,jsize))
  oscale(:,:)=0.
  do j =1,jsize
     do i =1,isize
        diffmax=-999.
        diffmin=999.
        do k =2,numpero*octaves
           if(diff(i,j,k).gt.diffmax) then
              diffmax=diff(i,j,k)
              dum_scale=k
           endif
           if(diff(i,j,k).lt.diffmin) then
              diffmin=diff(i,j,k)
           endif           
        enddo
!        if(confidence.le.0.0007) cycle
        if(diffmax.le.0.015) cycle
        oscale(i,j)=sigma*(k_int**real(dum_scale-1))*1.5
        confidence(i,j)=(diffmax+diffmin)*(diffmax-diffmin)
        if(oscale(i,j).le.2) cycle
        if(confidence(i,j).le.0.001) oscale(i,j)=0.
     enddo
  enddo
  !-----------ORIENTATION SELECT BY USING GRADIENT HISTGRAM---------------------
!  mask2(:,:)=0  
  allocate(orient(isize,jsize))
  call Orientation_Select(ndvi(:,:),orient(:,:),isize,jsize,oscale(:,:),mask2(:,:))

  !-----------CALCULATE SIMILARITIES BY A TEMPLATE MATCHING TECHNIQUE-------------  
   allocate(bbp(isize,jsize))

     CALL BBS(data(:,:,1:4),bbp(:,:),oscale(:,:),orient(:,:),mask2(:,:),isize,jsize,4,dir)
!!   mask2(:,:)=int(bbp(:,:)*1000)
!-----------------------WRITE THE RESULT DATA-------------------------------
   allocate(forshape(isize,jsize))
   forshape(:,:)=0
   do j =1,jsize
      do i =1,isize
         if(bbp(i,j).lt.0.4) cycle
            window=int(oscale(i,j)*2*1.5)
            mi=mod(window,2)
            if(mi.eq.0) window=window+1
            k=(window-1)/2
            print *,window
            conf_index=1
            do jj =1,window+4
               do ii =1,window+4
                  i2=i+ii-k-1-2
                  j2=j+jj-k-1-2
                  if(bbp(i2,j2).lt.0.4) cycle
                  if(confidence(i,j).lt.confidence(i2,j2)) then
                     conf_index=0;exit
                  endif
               enddo
               if(conf_index.eq.0) exit
            enddo

            if(conf_index.eq.1) then

               forshape(i-k:i+k,j-k:j+k)=1               
            endif
         enddo
      enddo
!  forshape(:,:)=data(:,:,1)
!      mask2(:,:)=forshape(:,:)
  CALL reshape_backforward_2D(forshape(:,:),data_o(:,:),isize,jsize)  
END SUBROUTINE LANDSLIDE_DETECTION


SUBROUTINE SAD(data_i,data_o,isize,jsize,imax_tl,jmax_tl,band)
    IMPLICIT NONE
  INTEGER :: i,j,k,ii,jj,kk,i2,j2,k2,i3,j3,k3,band  
  INTEGER(4),INTENT(IN) :: isize, jsize
  INTEGER(4),INTENT(IN) :: imax_tl,jmax_tl
  REAL(4),INTENT(OUT):: data_o(1:isize,1:jsize)
  REAL(4) :: templ(imax_tl,jmax_tl,band)
  REAL(4),INTENT(IN) :: data_i(isize,jsize,band)

!  CALL READ_TEMPLATE(templ,data_i,isize,jsize,imax_tl,jmax_tl,band)
  
  data_o(:,:)=0

!!$  do j =1+(jmax_tl-1)/2, jsize-(jmax_tl-1)/2
!!$     do i =1+(imax_tl-1)/2, isize-(imax_tl-1)/2
!!$        do jj=1,jmax_tl
!!$           do ii=1,imax_tl
!!$              i2=i+ii-1-(imax_tl-1)/2
!!$              j2=j+jj-1-(jmax_tl-1)/2
!!$           data_o(i,j)=abs(data_i(i2,j2)-templ(ii,jj))+data_o(i,j)
!!$       enddo
!!$       enddo
!!$  enddo
!!$  enddo
END SUBROUTINE SAD

  SUBROUTINE SSD(data_i,data_o,isize,jsize,imax_tl,jmax_tl,band)
    IMPLICIT NONE
  INTEGER :: i,j,k,ii,jj,kk,i2,j2,k2,i3,j3,k3    
  INTEGER(4),INTENT(IN) :: isize, jsize,band
  INTEGER(4),INTENT(IN) :: imax_tl,jmax_tl
  REAL(4),INTENT(OUT):: data_o(1:isize,1:jsize)
  REAL(4) :: templ(imax_tl,jmax_tl,band)
  REAL(4),INTENT(IN) :: data_i(isize,jsize,band)

!  CALL READ_TEMPLATE(templ,data_i,isize,jsize,imax_tl,jmax_tl,band)
!!$  do j =1,jmax_tl
!!$  do i =1,imax_tl
!!$     j2=j+int(real(j)*(-0.05))
!!$     if(j2.ge.1) then
!!$     templ(i,j)=templ(i,j2)
!!$!     print *,templ(i,j)
!!$     endif
!!$  enddo
!!$  enddo
  data_o(:,:)=0
!!$
!!$ 
!!$  do j =1+(jmax_tl-1)/2, jsize-(jmax_tl-1)/2
!!$     do i =1+(imax_tl-1)/2, isize-(imax_tl-1)/2
!!$        do jj=1,jmax_tl
!!$           do ii=1,imax_tl
!!$              i2=i+ii-1-(imax_tl-1)/2
!!$              j2=j+jj-1-(jmax_tl-1)/2
!!$           data_o(i,j)=(data_i(i2,j2,k)-templ(ii,jj,k))**2+data_o(i,j)
!!$       enddo
!!$       enddo
!!$  enddo
!!$  enddo
END SUBROUTINE SSD

  SUBROUTINE NCC(data_i,data_o,isize,jsize,imax_tl,jmax_tl,band)
    IMPLICIT NONE
  INTEGER :: i,j,k,ii,jj,kk,i2,j2,k2,i3,j3,k3,band
  INTEGER(4),INTENT(IN) :: isize, jsize
  INTEGER(4),INTENT(IN) :: imax_tl,jmax_tl
  REAL(4),INTENT(OUT):: data_o(1:isize,1:jsize)
  REAL(4) :: templ(imax_tl,jmax_tl)
  REAL(4),INTENT(IN) :: data_i(isize,jsize)
  REAL(4) :: xx,yy,xy
  
!  CALL READ_TEMPLATE(templ,data_i,isize,jsize,imax_tl,jmax_tl)
  
  data_o(:,:)=0

!!$  do j =1+(jmax_tl-1)/2, jsize-(jmax_tl-1)/2
!!$     do i =1+(imax_tl-1)/2, isize-(imax_tl-1)/2
!!$        xx=0.;yy=0.;xy=0.
!!$        do jj=1,jmax_tl
!!$           do ii=1,imax_tl
!!$              i2=i+ii-1-(imax_tl-1)/2
!!$              j2=j+jj-1-(jmax_tl-1)/2
!!$              xx=data_i(i2,j2)**2+xx
!!$              yy=templ(ii,jj)**2+yy
!!$              xy=data_i(i2,j2)*templ(ii,jj)+xy
!!$       enddo
!!$    enddo
!!$    data_o(i,j)=xy/sqrt(xx*yy)
!!$  enddo
!!$  enddo
END SUBROUTINE NCC

  SUBROUTINE ZNCC(data_i,data_o,oscale,orient,isize,jsize,band,dir)
    IMPLICIT NONE
    INTEGER :: i,j,ii,jj,k,kk,iii,jjj,i2,j2,i3,j3 
  REAL(8) :: x1,x2,y1,y2,xx,yy,xy
  INTEGER(4) :: inum
  INTEGER(4),INTENT(IN) :: isize, jsize,band
  INTEGER(4),PARAMETER :: imax_tl=17
  INTEGER(4),PARAMETER :: jmax_tl=17
  REAL(4),INTENT(OUT):: data_o(1:isize,1:jsize)  
  REAL(4) :: templ(imax_tl,jmax_tl,band)
  REAL(4),INTENT(IN) :: data_i(isize,jsize,band)
  REAL(4),INTENT(IN) :: oscale(isize,jsize)
  REAL(4),INTENT(IN) :: orient(isize,jsize)
  INTEGER(4) :: ii1,ii2,jj1,jj2,kx,ky,window
  INTEGER(4) :: mi,mj,isize2,jsize2,rad
  REAL(4),ALLOCATABLE :: xx1(:,:),xx2(:,:),xx3(:,:),xx4(:,:)
  REAL(4),ALLOCATABLE :: yy1(:,:),yy2(:,:),yy3(:,:),yy4(:,:)
  REAL(4) :: theta
  character(512),INTENT(IN) :: dir
  CALL READ_TEMPLATE(templ,imax_tl,jmax_tl,band,2,dir)
!  print *,'hi'  
  kx=(imax_tl-1)/2
  ky=(jmax_tl-1)/2
  theta=0
  allocate(xx3(imax_tl,jmax_tl),xx4(imax_tl,jmax_tl))
  allocate(yy3(imax_tl,jmax_tl),yy4(imax_tl,jmax_tl))
 
        do jj = 1, jmax_tl
           do ii = 1, imax_tl
!              print *,ii,jj
              i2=ii-kx-1
              j2=jj-ky-1
           xx3(ii,jj)=real(i2)*cos(theta)+real(j2)*sin(theta)
           yy3(ii,jj)=real(j2)*cos(theta)-real(i2)*sin(theta)
           xx4(ii,jj)=xx3(ii,jj)/real(imax_tl-1)
           yy4(ii,jj)=yy3(ii,jj)/real(jmax_tl-1)
!           print *,x3(ii,jj),y3(ii,jj)
        enddo
       enddo

  data_o(:,:)=0

  do j =1+(jmax_tl-1)/2, jsize-(jmax_tl-1)/2
     do i =1+(imax_tl-1)/2, isize-(imax_tl-1)/2  
        if(oscale(i,j).eq.0) cycle
!        print *,i,j
!-----------------initialization---------------------------
        window=nint(oscale(i,j)*2*1.5)    !1.5 times radius 
        mi=mod(window,2)
        if(mi.eq.0) window=window+1
        rad=(window-1)/2
        theta=0. !orient(i,j)*dora
!        print *,window,k
        allocate(xx1(window,window),yy1(window,window))
        allocate(xx2(window,window),yy2(window,window))
        do jj = 1, window
           do ii = 1, window
              i2=ii-rad-1
              j2=jj-rad-1
              xx1(ii,jj)=real(i2)*cos(theta)+real(j2)*sin(theta)
              yy1(ii,jj)=real(j2)*cos(theta)-real(i2)*sin(theta)
              xx2(ii,jj)=xx1(ii,jj)/real(window-1)
              yy2(ii,jj)=yy1(ii,jj)/real(window-1)
!              print *,x1(ii,jj),x2(ii,jj)
           enddo
        enddo
        !----------------------estimate CCs-------------------------------
        do k =1,band
        x1=0.;x2=0.;y1=0.;y2=0.;xy=0.
        inum=0
        do jj=1,jmax_tl
           do ii=1,imax_tl
              inum=inum+1
              i2=i+ii-1-(imax_tl-1)/2
              j2=j+jj-1-(jmax_tl-1)/2
              x1=data_i(i2,j2,k)+x1
              x2=data_i(i2,j2,k)**2+x2
              y1=templ(ii,jj,k)+y1
              y2=templ(ii,jj,k)**2+y2              
              xy=data_i(i2,j2,k)*templ(ii,jj,k)+xy
       enddo
    enddo
!    print *,x1,x2,y1,y2,xy
    xx=real(x2)/real(inum)-(real(x1)/real(inum))**2
    yy=real(y2)/real(inum)-(real(y1)/real(inum))**2
    xy=real(xy)/real(inum)-real(x1*y1)/real(inum**2)
    data_o(i,j)=xy/sqrt(xx*yy)+data_o(i,j)
 enddo

        deallocate(xx1,yy1,xx2,yy2)
 
      enddo
  enddo
END SUBROUTINE ZNCC

  SUBROUTINE BBS(data_i,bbp,oscale,orient,mask,isize,jsize,band,dir)
    IMPLICIT NONE
  INTEGER :: i,j,ii,jj,k,kk,iii,jjj,i2,j2,i3,j3
  INTEGER(4),INTENT(IN) :: isize, jsize,band
  INTEGER(4),PARAMETER :: imax_tl=17
  INTEGER(4),PARAMETER :: jmax_tl=17
  INTEGER(4),PARAMETER :: n_templ=2
  REAL(4),INTENT(OUT):: bbp(1:isize,1:jsize)
  REAL(4) :: templ(imax_tl,jmax_tl,band,n_templ)
  REAL(4),INTENT(IN) :: data_i(isize,jsize,band)
  REAL(4),INTENT(IN) :: oscale(isize,jsize)
  REAL(4),INTENT(IN) :: orient(isize,jsize)
  INTEGER(4),INTENT(IN) :: mask(isize,jsize)
  REAL(4),ALLOCATABLE :: x1(:,:),x2(:,:),x3(:,:,:),x4(:,:,:)
  REAL(4),ALLOCATABLE :: y1(:,:),y2(:,:),y3(:,:,:),y4(:,:,:)
  REAL(4) :: ai1,ai2,aj1,aj2,dis_rgb,dis_xy,pi,dora,theta
  INTEGER(4) :: ii1,ii2,jj1,jj2,kx,ky,window
  INTEGER(4) :: inum
  REAL(4),parameter :: disdep= 0.01  !dependency factor of pixels distance
  INTEGER(4),PARAMETER :: psize=1
  INTEGER(4),ALLOCATABLE :: minx(:,:), miny(:,:)
  INTEGER(4) :: minx2(imax_tl,jmax_tl), miny2(imax_tl,jmax_tl)
  REAL(4) :: distance,maxdis,adis,bbp_max,asphalt
  INTEGER(4) :: mi,mj,isize2,jsize2,rad,nt
  real(4) :: dis(imax_tl,jmax_tl)
  character(512),INTENT(IN) :: dir
  
  pi=acos(-1.)
  dora=pi/180.
  bbp(:,:)=0

  CALL READ_TEMPLATE(templ,imax_tl,jmax_tl,band,n_templ,dir)
  kx=(imax_tl-1)/2
  ky=(jmax_tl-1)/2
  theta=0
  allocate(x3(imax_tl,jmax_tl,n_templ),x4(imax_tl,jmax_tl,n_templ))
  allocate(y3(imax_tl,jmax_tl,n_templ),y4(imax_tl,jmax_tl,n_templ))
  do nt=1,n_templ
  maxdis=sqrt(real(kx)**2+real(ky)**2)  
        do jj = 1, jmax_tl
           do ii = 1, imax_tl
              i2=ii-kx-1
              j2=jj-ky-1
           x3(ii,jj,nt)=real(i2)*cos(theta)+real(j2)*sin(theta)
           y3(ii,jj,nt)=real(j2)*cos(theta)-real(i2)*sin(theta)
           x4(ii,jj,nt)=x3(ii,jj,nt)/real(imax_tl-1)
           y4(ii,jj,nt)=y3(ii,jj,nt)/real(jmax_tl-1)
!           print *,x3(ii,jj),y3(ii,jj)
        enddo
       enddo
   enddo

  do j =1, jsize
     do i =1, isize
!        print *,i,j
        if(mask(i,j).eq.999.) cycle
        if(oscale(i,j).le.1.5) cycle
!-----------------initialization---------------------------
        window=nint(oscale(i,j)*2.0*2.0)    !1.5 times radius 
        mi=mod(window,2)
        if(mi.eq.0) window=window+1
        rad=(window-1)/2
!        theta=(orient(i,j)-145.)*dora
!        print *,window,k
        allocate(x1(window,window),y1(window,window))
        allocate(x2(window,window),y2(window,window))
        allocate(minx(window,window),miny(window,window))


        bbp_max=0
        do nt =1, n_templ
        do kk =1,2        
           if(kk.eq.1) theta=(orient(i,j)-145.)*dora
           if(kk.eq.2) theta=((orient(i,j)-145.)-180)*dora
        
        do jj = 1, window
           do ii = 1, window
              i2=ii-rad-1
              j2=jj-rad-1
              x1(ii,jj)=real(i2)*cos(theta)+real(j2)*sin(theta)
              y1(ii,jj)=real(j2)*cos(theta)-real(i2)*sin(theta)
              x2(ii,jj)=x1(ii,jj)/real(window-1)
              y2(ii,jj)=y1(ii,jj)/real(window-1)
!              print *,x1(ii,jj),x2(ii,jj)
           enddo
        enddo
!----------------------estimate BBPs-------------------------------
        inum=0
        minx(:,:)=999
        miny(:,:)=999
        dis(:,:)=10000.

        do jj=1,window
           do ii=1,window
              distance=1000000.
              ii1=int(x1(ii,jj))+i;jj1=int(y1(ii,jj))+j
              ai1=x2(ii,jj);aj1=y2(ii,jj)
              if((ii1.lt.1).or.(ii1.gt.isize).or.(jj1.lt.1).or.(jj1.gt.jsize)) cycle
              
              do jjj=1,jmax_tl,psize
                 do iii=1,imax_tl,psize
                    ii2=int(x3(iii,jjj,nt));jj2=int(y3(iii,jjj,nt))
                    ai2=x4(iii,jjj,nt);aj2=y4(iii,jjj,nt)

                    dis_xy=real(ai1-ai2)**2+real(aj1-aj2)**2
                    dis_rgb=0
                    do k=1,band
                       dis_rgb=(data_i(ii1,jj1,k)/65536.-templ(iii,jjj,k,nt)/65536.)**2 + dis_rgb
                    enddo
                    adis=dis_rgb+dis_xy*disdep
 !                   print *,dis_rgb,dis_xy*disdep
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
              enddo

           enddo
        enddo
 
        do jj=1,window
           do ii=1,window
              if(minx(ii,jj).eq.999) cycle
              i2=minx(ii,jj)
              j2=miny(ii,jj)
              i3=minx2(i2,j2)
              j3=miny2(i2,j2)

              ii1=int(x1(i3,j3))+i;jj1=int(y1(i3,j3))+j
              asphalt=data_i(ii1,jj1,3)-(data_i(ii1,jj1,1)*2.-11200)
              if(asphalt.le.0) cycle
              if((ii.eq.i3).and.(jj.eq.j3)) then
                 bbp(i,j)=bbp(i,j)+1
              endif
              
           enddo
        enddo
        
           bbp(i,j)=bbp(i,j)/min(window*window,imax_tl*jmax_tl/(psize**2))
        if(bbp(i,j).ge.bbp_max) then
           bbp_max=bbp(i,j)
        endif
        
     enddo
     enddo

     bbp(i,j)=bbp_max
     deallocate(x1,y1,x2,y2,minx,miny)

     
     enddo
  enddo
END SUBROUTINE BBS

  SUBROUTINE Orientation_Select(data_i,data_o,isize,jsize,oscale,mask)
    IMPLICIT NONE
  INTEGER :: i,j,k,ii,jj,kk,i2,j2,k2,i3,j3,k3    
  INTEGER(4),INTENT(IN) :: isize, jsize
  REAL(4),INTENT(OUT):: data_o(1:isize,1:jsize)
  REAL(4),INTENT(IN) :: data_i(isize,jsize)
  INTEGER(4),INTENT(IN) :: mask(isize,jsize)
  REAL(4),INTENT(IN) :: oscale(isize,jsize)
  INTEGER(4),PARAMETER :: interval = 10
  REAL(4) :: hist(180/interval)
  REAL(4) :: dx,dy,pi,dora,theta,hismax,mu
  INTEGER :: itheta,window,mi

  pi=acos(-1.)
  dora=pi/180
  data_o(:,:)=999.

  do j =1, jsize
     do i =1, isize
        if(mask(i,j).eq.999) cycle
        window=nint(oscale(i,j)*2.0*1.5)    !1.5 times radius 
        mi=mod(window,2)
        if(mi.eq.0) window=window+1
        k=(window-1)/2
        if((window.lt.3)) cycle
        if((i+k+1.gt.isize).or.(i-k-1.lt.1)) cycle
        if((j+k+1.gt.jsize).or.(j-k-1.lt.1)) cycle
        hist(:)=0
        do jj=1,window
           do ii=1,window
              i2=i+ii-1-k
              j2=j+jj-1-k
              dx=data_i(i2+1,j2)-data_i(i2-1,j2)
              dy=data_i(i2,j2+1)-data_i(i2,j2-1)
              if(dx.eq.0) cycle
              mu=dx**2+dy**2
              theta=atan(dy/dx)/dora
!              print *,theta,dy,dx
              if(theta.le.0) theta=theta+180
              itheta=int(theta/interval)+1
              hist(itheta)=mu+hist(itheta)
           enddo
        enddo

        hismax=0
        do ii = 1,180/interval
           hismax=max(hist(ii),hismax)
           if(hismax.eq.hist(ii)) then
              mu=ii*10-5
              if(mu.le.90) data_o(i,j)=180.-(mu+90)
              if(mu.ge.90) data_o(i,j)=180.-(mu-90)
           endif
       enddo
        
     enddo
  enddo
END SUBROUTINE ORIENTATION_SELECT

  SUBROUTINE READ_TEMPLATE(tmpl,isize2,jsize2,band2,n_templ,dir)
    IMPLICIT NONE
  INTEGER :: i,j,k,isize,jsize,band
  INTEGER(4),INTENT(IN) :: isize2,jsize2,band2,n_templ
  REAL(4),INTENT(OUT) :: tmpl(isize2,jsize2,band2,n_templ)
  character(512) :: infile
  character(512),INTENT(IN) :: dir
  infile='template.tif'
  infile=trim(dir)//infile
  call read_size_tif(infile,isize,jsize,band)
  CALL open_tif(infile,tmpl(:,:,:,1),isize,jsize,band,jsize)
  if(n_templ.eq.2) then
  infile='template2.tif'
  infile=trim(dir)//infile
  call read_size_tif(infile,isize,jsize,band)  
  CALL open_tif(infile,tmpl(:,:,:,2),isize,jsize,band,jsize)
  endif
!  tmpl(:,:,:,:)=1
END SUBROUTINE READ_TEMPLATE

END MODULE Template_Matching_lr
