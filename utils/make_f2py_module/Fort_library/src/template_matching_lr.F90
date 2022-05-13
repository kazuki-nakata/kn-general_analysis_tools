MODULE Template_Matching_lr
  USE Basic_processing
  IMPLICIT NONE
  INTEGER, PARAMETER :: MAX_IO = 64

CONTAINS

  SUBROUTINE landslide_detection(data_i,data_o,band,jsize,isize,tmp_i,num_t,band_t,jsize_t,isize_t,ltfile&
       ,ndvi_thres,bbs_thres,confidence_thres,scale_thres)
  
  INTEGER :: i,j,k,l,ii,jj,kk,i2,j2,k2,i3,j3,k3      
  REAL, PARAMETER :: threshold_b= 0.4
  INTEGER, PARAMETER :: band1=3
  INTEGER, PARAMETER :: band2=4
  INTEGER(4),INTENT(IN) :: isize,jsize,band
  INTEGER(4),INTENT(IN) :: isize_t,jsize_t,band_t,num_t
  REAL(4),INTENT(IN) :: ndvi_thres,bbs_thres,confidence_thres,scale_thres
  INTEGER(4), INTENT(OUT)  :: data_o(jsize,isize)
  REAL(4), INTENT(IN) :: data_i(band,jsize,isize)
  REAL(4), INTENT(IN) :: tmp_i(num_t,band_t,jsize_t,isize_t)
  INTEGER(4), ALLOCATABLE :: morp1(:,:),morp2(:,:),bin(:,:)
  REAL(4), ALLOCATABLE :: dum(:,:),ndvi(:,:),data(:,:,:)
  INTEGER(4), ALLOCATABLE :: mask1(:,:),mask2(:,:),mask3(:,:)
  REAL(4), ALLOCATABLE :: octave(:,:,:),oscale(:,:),orient(:,:)
  REAL(4), ALLOCATABLE :: diff(:,:,:),diff2(:,:,:)
  REAL(4), ALLOCATABLE :: dummy(:,:),confidence(:,:),bbp(:,:)
  INTEGER(4), ALLOCATABLE :: dummy2(:,:),forshape(:,:)
  INTEGER(1), ALLOCATABLE :: omax_scale(:,:,:)
  REAL(4), ALLOCATABLE :: tmp0(:,:,:,:)
  REAL(4), ALLOCATABLE :: tmp_orient(:)
  INTEGER(4), PARAMETER :: imax_tl=17 !17
  INTEGER(4), PARAMETER :: jmax_tl=17 !17
  INTEGER(4) :: idirect,jdirect
  REAL(4), ALLOCATABLE :: tmp(:,:,:,:)
  !------test data---
  REAL(4), ALLOCATABLE :: test(:,:,:)
  !---------
  INTEGER, PARAMETER :: octaves = 2
  INTEGER, PARAMETER :: interval = 2
  INTEGER, PARAMETER :: numpero = 6
  REAL, PARAMETER :: sigma=0.5
  REAL(4) :: sum,ksigma,k_int,diffmax,sf,diffmin,diff_dum,sum_conf,sum_conf2
  INTEGER :: mi,mj,wsize,dum_scale,disize,djsize,n_oc,interval2,first,sum_max
  INTEGER(4) :: isize2,jsize2,isize3,jsize3,conf_index,window,diffmax_scale
  REAL(4) ::pi,dora,theta
  REAL(4) :: gain(4),offset(4),abscalfactor(4),effectivebandwidth(4)
  REAL(4) :: cal_gain,cal_offset
  character(4) :: ti
  character(512),INTENT(IN) :: ltfile
  character(512) ::infile,dir
  REAL(4) ::dist,overw

  pi=acos(-1.)
  dora=pi/180.


  print *,'hello'
  print *, bbs_thres,ndvi_thres
  !---------------INITIALIZATION-----------------------------
  allocate(ndvi(isize,jsize),bin(isize,jsize))  
  allocate(mask1(isize,jsize),mask2(isize,jsize),mask3(isize,jsize))
  allocate(data(isize,jsize,band))
  allocate(test(isize,jsize,3))    
  allocate(tmp0(isize_t,jsize_t,band_t,num_t))
  allocate(tmp(imax_tl,jmax_tl,band_t,num_t))
  allocate(tmp_orient(num_t))
  CALL reshape_backforward_4D(tmp_i,tmp0,num_t,band_t,jsize_t,isize_t)
  CALL ARRANGE_TEMPLATE(tmp,imax_tl,jmax_tl,band_t,num_t,tmp0,isize_t,jsize_t,tmp_orient)
  CALL reshape_backforward_3D(data_i(1:band,1:jsize,1:isize),data(1:isize,1:jsize,1:band),band,jsize,isize)

!!$  
!!$  allocate(test(isize,jsize,band))
!!$  test(:,:,:)=data(:,:,:)
!!$  data(:,:,3)=test(:,:,1)
!!$  data(:,:,2)=test(:,:,2)
!!$  data(:,:,1)=test(:,:,3)

  do j =1,jsize
     do i=1,isize
        if(data(i,j,2).eq.0) then
           mask3(i,j)=999
        endif
     enddo
  enddo
  
  !mask1=dynamic pyramid nodata edge mask, mask2= all mask data, mask3=all mask data
!  CALL reshape_backforward_2D(tmp0(:,:,2,1),data_o(:,:),isize,jsize)

  !--------------MAKE MASK DATA USING LAND COVER DATA-----------------------
  CALL Index(data(1:isize,1:jsize,1:band),ndvi(1:isize,1:jsize),isize,jsize,band,band1,band2,0)
  CALL Binarization(ndvi(1:isize,1:jsize),bin(1:isize,1:jsize),isize,jsize,ndvi_thres,999)
  CALL Make_Mask_FORLANDSLIDE(data(:,:,4:5),mask1(:,:),isize,jsize)
  CALL Mask_MERGE(bin(:,:),mask1(:,:),mask2(:,:),isize,jsize)
  deallocate(bin)

  CALL Make_MASK_NODATA(ndvi(:,:),mask3(:,:),isize,jsize,999)
  mask1(:,:)=mask3(:,:)

  do j =1,jsize
     do i=1,isize
        if((data(i,j,1).ge.0.99).and.(data(i,j,2).ge.0.99).and.(data(i,j,3).ge.0.99)) then
           mask2(i,j)=999;mask3(i,j)=999
        endif
     enddo
  enddo
  
!!$  dir='D:\\JR\\mask.tif'
!!$  test(:,:,1)=mask1(:,:)
!!$  test(:,:,2)=mask2(:,:)
!!$  test(:,:,3)=mask3(:,:)
!!$  call write_tif(dir,test,isize,jsize,3,jsize)  

  !--------------MAKE MASK DATA FOR SMALL SCALE REGION----------
  !---opening for removing small scale objects in high ndvi regions
  !---closing for removing small scale objects in low ndvi regions
   
!!$  CALL Index(data(:,:,:),ndvi(:,:),isize,jsize,band,band1,band2,0)
!!$  CALL Interpolation_type1(ndvi(:,:),mask1(:,:),isize,jsize,999.)
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

!!$  !---------------CREATE SCALE-SPACE MAPS BY DOG OPERATOR----------  

  wsize=6*interval+1
  k_int=(real(interval)/sigma)**(1/real(numpero))
  allocate(diff(isize,jsize,octaves),diff2(isize,jsize,octaves))
  allocate(omax_scale(isize,jsize,octaves))  
  diff(:,:,:)=0.
  diff2(:,:,:)=0.
  omax_scale(:,:,:)=0
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

     !     print *, real(interval2)/sf

     do i =1, numpero
        if(n_oc.eq.1) ksigma = sigma*(k_int**real(i))
        if(n_oc.ne.1) ksigma= sqrt((real(interval2)*(k_int**real(i)))**2-real(interval2)**2)
        print *,ksigma,sqrt(ksigma**2+interval2**2)
        CALL gaussian_filter(octave(:,:,0),octave(:,:,i),mask1(:,:),isize2,jsize2,wsize,ksigma,sf) 
     enddo

     if(n_oc.eq.1) first=2
     if(n_oc.eq.2) first=1
     kk=n_oc   
     do j =1,jsize
        do i =1,isize
           if(mask3(i,j).eq.999) cycle
           diffmax=-999.
           diffmin=999.
           do k =first,numpero
              i2=(i-1)/interval2+1
              j2=(j-1)/interval2+1
              if((i2.gt.isize2).or.(j2.gt.jsize2)) cycle
              diff_dum=octave(i2,j2,k)-octave(i2,j2,k-1)
              if(diff_dum.ge.diffmax) then
                 diffmax=diff_dum
                 diffmax_scale=(kk-1)*numpero+k
              endif
              if(diff_dum.le.diffmin) diffmin=diff_dum          
           enddo
           diff(i,j,kk)=diffmax
           omax_scale(i,j,kk)=diffmax_scale
           diff2(i,j,kk)=diffmin
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


  deallocate(dummy,dummy2)
  
  print *,'hi'
  !------scale selection by searching a maximum value----------
  allocate(confidence(isize,jsize))
  confidence(:,:)=0.
  allocate(oscale(isize,jsize))
  oscale(:,:)=0.

  do j =1,jsize
     do i =1,isize
        diffmax=-999.
        diffmin=999.
        do k =1,octaves
           if(diff(i,j,k).gt.diffmax) then
              diffmax=diff(i,j,k)
              dum_scale=int(omax_scale(i,j,k))
!              print *,dum_scale
 !          print *,diff(i,j,k),omax_scale(i,j,k),i,j
           endif
           if(diff2(i,j,k).lt.diffmin) then
              diffmin=diff2(i,j,k)
           endif           
        enddo
        if(diffmax.le.0.01) cycle
        oscale(i,j)=sigma*(k_int**real(dum_scale-1))*1.5
        confidence(i,j)=(diffmax+diffmin)*(diffmax-diffmin)
     enddo
  enddo
 !  print *,'hi'

!!$  dir='D:\\JR\\diff.tif'
!!$  call write_tif(dir,diff,isize,jsize,7,jsize)
  
  deallocate(diff,diff2,omax_scale)

  !-------------select (mask) candidate windows------------------
  !-----1. confidence mask
  !-----2. mask unnecessary window using scale (overwrap) and confidence
  !-----3. scale mask (this mask must be put after the above mask 2)
  
  do j =1, jsize
     do i =1, isize
!        print *,mask2(i,j),confidence(i,j),oscale(i,j)
        if(mask2(i,j).eq.999) cycle
        if(confidence(i,j).le.0.0002) then     
           mask2(i,j)=999
           oscale(i,j)=0.
        endif
     enddo
  enddo
  
  do l =1,2
     if(l.eq.1) then
        idirect=1
        jdirect=1
     elseif(l.eq.2) then
        idirect=-1
        jdirect=1
     endif
     
     do j =1, jsize,jdirect
     do i =1, isize,idirect
        if(mask2(i,j).eq.999) cycle

        window=nint(oscale(i,j)*2)
        mi=mod(window,2)
        if(mi.eq.0) window=window+1
        k=(window-1)/2  

        if((i+k.gt.isize).or.(i-k.lt.1)) cycle
        if((j+k.gt.jsize).or.(j-k.lt.1)) cycle

        do jj=1,window
           do ii=1,window
              i2=i+ii-1-k
              j2=j+jj-1-k
              if((i.eq.i2).and.(j.eq.j2)) cycle
              
              if(l.eq.1) then
                 if(mask2(i2,j2).eq.999) cycle
              endif
              
              dist=sqrt(real((i-i2)**2+(j-j2)**2))
              overw=real(oscale(i,j)/2.+oscale(i2,j2)/2.-dist)/real(oscale(i,j)/2.)
              !              print *,overw
              if(overw.lt.0.66) cycle
              if(confidence(i,j).lt.confidence(i2,j2)) then
                 mask2(i,j)=999;exit
              else
                 mask2(i2,j2)=999
              endif                
           enddo
        enddo
     enddo
  enddo
enddo

  do j =1, jsize
     do i =1, isize
        if(mask2(i,j).eq.999) cycle
         if(oscale(i,j).lt.scale_thres) then      
           mask2(i,j)=999
        endif
     enddo
  enddo
  
!!$  test(:,:,1)=confidence(:,:)
!!$  test(:,:,2)=oscale(:,:)
!!$  test(:,:,3)=data(:,:,1)
!!$  dir='D:\\JR\\conf.tif'
!!$  call write_tif(dir,test,isize,jsize,3,jsize)


  !------------------------test--------------------------------
!!$  allocate(orient(isize,jsize))
!!$  call Orientation_Select(ndvi(:,:),orient(:,:),isize,jsize,oscale(:,:),mask2(:,:))
!!$   allocate(forshape(isize,jsize))
!!$   forshape(:,:)=0
!!$   
!!$   do j =1,jsize
!!$      do i =1,isize
!!$         if(mask2(i,j).eq.999) cycle
!!$            window=int(oscale(i,j))
!!$            mi=mod(window,2)
!!$            if(mi.eq.0) window=window+1
!!$            k=(window-1)/2
!!$        if((i+k.gt.isize).or.(i-k.lt.1)) cycle
!!$        if((j+k.gt.jsize).or.(j-k.lt.1)) cycle
!!$               forshape(i-k:i+k,j-k:j+k)=1
!!$         enddo
!!$      enddo
!!$      
!!$  print *,'writing'
!!$  test(:,:,1)=real(data(:,:,4))
!!$  test(:,:,2)=real(forshape(:,:))
!!$  test(:,:,3)=real(orient(:,:))
!!$  dir='D:\\JR\\mask3.tif'
!!$  call write_tif(dir,test,isize,jsize,3,jsize)
!!$  print *,'finish'
!!$  
!!$  CALL reshape_backforward_2D(forshape(:,:),data_o(:,:),isize,jsize)

  
!---------------------------------------------------------------
  !-----------ORIENTATION SELECT BY USING GRADIENT HISTGRAM---------------------
  allocate(orient(isize,jsize))
  call Orientation_Select(ndvi(:,:),orient(:,:),isize,jsize,oscale(:,:),mask2(:,:))
  
  print *,'orientation finish'
  !-----------CALCULATE SIMILARITIES BY A TEMPLATE MATCHING TECHNIQUE-------------  
   allocate(bbp(isize,jsize))
   CALL BBS(data(:,:,1:4),bbp(:,:),oscale(:,:),orient(:,:),mask2(:,:),isize,jsize,4, &
        tmp(1:imax_tl,1:jmax_tl,1:4,:),imax_tl,jmax_tl,num_t,tmp_orient,ltfile,bbs_thres)
!!   mask2(:,:)=int(bbp(:,:)*1000)
   !-----------------------WRITE THE RESULT DATA-------------------------------
!!$   allocate(forshape(isize,jsize))
!!$   forshape(:,:)=0
!!$   
!!$   do j =1,jsize
!!$      do i =1,isize
!!$         if(mask2(i,j).eq.999) cycle
!!$         if(bbp(i,j).lt.0.3) cycle
!!$            window=int(oscale(i,j)*2*1.5)
!!$            mi=mod(window,2)
!!$            if(mi.eq.0) window=window+1
!!$            k=(window-1)/2
!!$        if((i+k.gt.isize).or.(i-k.lt.1)) cycle
!!$        if((j+k.gt.jsize).or.(j-k.lt.1)) cycle
!!$               forshape(i-k:i+k,j-k:j+k)=1
!!$         enddo
!!$      enddo
!------------------------------------------------------------------------------
   allocate(forshape(isize,jsize))
   forshape(:,:)=0
   
   do j =1,jsize
      do i =1,isize
         if(mask2(i,j).eq.999) cycle
         if(bbp(i,j).lt.0.0) cycle
            window=int(oscale(i,j)*2)
            mi=mod(window,2)
            if(mi.eq.0) window=window+1
            k=(window-1)/2
            conf_index=1
            do jj =1,window
               do ii =1,window
                  i2=i+ii-k-1-2
                  j2=j+jj-k-1-2
                  if(bbp(i2,j2).lt.0.0) cycle
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
!------------------test------------------------------      
!!$  print *,'writing'
!!$  test(:,:,1)=real(bbp(:,:))
!!$  test(:,:,2)=real(forshape(:,:))
!!$  test(:,:,3)=real(orient(:,:))
!!$  dir='D:\\JR\\mask3.tif'
!!$  call write_tif(dir,test,isize,jsize,3,jsize)
!!$  print *,'finish'
  
  CALL reshape_backforward_2D(forshape(:,:),data_o(:,:),isize,jsize)

END SUBROUTINE landslide_detection

  SUBROUTINE BBS(data_i,bbp,oscale,orient,mask,isize,jsize,band,templ,imax_tl,jmax_tl,n_templ,tmp_orient,ltfile,bbs_thres)
    IMPLICIT NONE
  INTEGER :: i,j,ii,jj,k,kk,iii,jjj,i2,j2,i3,j3,i4,j4,scan_i,scan_j
  INTEGER(4) :: ii1,ii2,jj1,jj2,kx,ky,window
  INTEGER(4) :: inum
  INTEGER(4) :: mi,mj,isize2,jsize2,rad,nt,step      
  INTEGER(4),INTENT(IN) :: isize, jsize,band
  INTEGER(4),INTENT(IN) :: imax_tl,jmax_tl,n_templ
  REAL(4) :: templ(imax_tl,jmax_tl,band,n_templ)
  INTEGER(4) :: minx2(imax_tl,jmax_tl), miny2(imax_tl,jmax_tl)  
  REAL(4),INTENT(OUT):: bbp(1:isize,1:jsize)
  REAL(4),INTENT(IN) :: data_i(isize,jsize,band)
  REAL(4),INTENT(IN) :: oscale(isize,jsize)
  REAL(4),INTENT(IN) :: orient(isize,jsize)
  REAL(4),INTENT(IN) :: tmp_orient(n_templ)
  INTEGER(4),INTENT(IN) :: mask(isize,jsize)
  INTEGER(4),ALLOCATABLE :: minx(:,:), miny(:,:)  
  REAL(4),ALLOCATABLE :: x1(:,:),x2(:,:),x3(:,:,:),x4(:,:,:)
  REAL(4),ALLOCATABLE :: y1(:,:),y2(:,:),y3(:,:,:),y4(:,:,:)
  real(4),ALLOCATABLE :: ndvi1(:,:,:),si1(:,:,:),conv(:,:,:)
  REAL(8),PARAMETER :: disdep= 4!0.2  !dependency factor of pixels distance
  INTEGER(4),PARAMETER :: psize=1
  REAL(4) :: ai1,ai2,aj1,aj2,pi,dora,theta
  REAL(8) :: distance,maxdis,adis,dis_rgb,dis_xy
  real(4) :: bbp_max,asphalt
  real(8) :: dis(imax_tl,jmax_tl)
  real(4) :: sum,arad
  integer(4) :: temp_sum  
  real(4) :: ndvi2,si2
  !-------------for normalization----------------
  real(4) :: ave(4,n_templ),ave2(4)
  real(4) :: a,b
  INTEGER(4) :: inum2,imax,imin
  real(4),ALLOCATABLE :: list(:)
  !--------------threshold--------------------
  real(4),intent(in) :: bbs_thres
  integer(4) :: isc,istd
  real(4) :: bbref
  character(512),INTENT(IN) :: ltfile
  real(4) :: bbs_ref0(40,10)
  real(4) :: bbs_ref00(800,10)
  real(4) :: bbs_ref(800)
  real(4) :: bbs_ref_min
  real(4) :: a1,a2,sum_min
  print *,ltfile
  open(70,file=ltfile,access='sequential',status='old',form='formatted')
  do j=1,40
  do i=1,10
     read(70,*) isc,istd,bbref
     bbs_ref0(j,i)=bbref
     print *,isc,istd,bbref,i,j
  enddo
  enddo
  do j=1,10
     do i=4,800
        ii1=(i-4)/10+1
        ii2=(i-4)/10+2
        a1=real(ii1-1)*10+4
        a2=real(ii2-1)*10+4
        a=(real(i)-a1)/(a2-a1)
        if(ii1.lt.1) ii1=1
        if(ii2.lt.1) ii2=1
        if(ii1.gt.40) ii1=40
        if(ii2.gt.40) ii2=40
        bbs_ref00(i,j)=(1-a)*bbs_ref0(ii1,j)+a*bbs_ref0(ii2,j)
     enddo
  enddo
  
  bbs_ref_min=10000.
  do ii=1,10
     if(abs(bbs_thres-bbs_ref00(289,ii)).le.bbs_ref_min) then
        bbs_ref_min=abs(bbs_thres-bbs_ref00(289,ii))
        jj=ii
     endif
  enddo
  bbs_ref(:)=bbs_ref00(:,jj)
!!$  do i=1,800
!!$     print *, bbs_ref(i),i,jj
!!$  enddo
  
  pi=acos(-1.)
  dora=pi/180.
  bbp(:,:)=-1
  
  kx=(imax_tl-1)/2
  ky=(jmax_tl-1)/2
  theta=0
  allocate(x3(imax_tl,jmax_tl,n_templ),x4(imax_tl,jmax_tl,n_templ))
  allocate(y3(imax_tl,jmax_tl,n_templ),y4(imax_tl,jmax_tl,n_templ))
  allocate(ndvi1(imax_tl,jmax_tl,n_templ),si1(imax_tl,jmax_tl,n_templ))

!-----------calc distance (-0.5 - 0.5). must set theta=0 (;template image projection base)------
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
        enddo
     enddo
  enddo
  
  temp_sum=0     
    do jj = 1, jmax_tl
       do ii = 1, imax_tl
          if((x4(ii,jj,1).lt.-0.51).or.(x4(ii,jj,1).gt.0.51).or. &
               (y4(ii,jj,1).lt.-0.51).or.(y4(ii,jj,1).gt.0.51)) cycle
          temp_sum=temp_sum+1
       enddo
    enddo
    
    print *,"test1=",imax_tl*jmax_tl,temp_sum  
    allocate(list(temp_sum))
    do nt=1,n_templ
       do k =1,4
          inum2=0          
          do jj = 1, jmax_tl
             do ii = 1, imax_tl
                inum2=1+inum2
                list(inum2)=templ(ii,jj,k,nt)
             enddo
          enddo
          call heapsort2(inum2,list)
          imax=nint(0.97*real(inum2))
          imin=nint(0.03*real(inum2))
          a=1/(list(imax)-list(imin))
          b=list(imin)
!!$          print *,a,b,k
!!$          print *,templ(:,10,k,nt)
          templ(:,:,k,nt)=a*(templ(:,:,k,nt)-b)
!          print *,templ(:,10,k,nt)
       enddo
    enddo
    deallocate(list)
    
!---------------------- main loop -----------------------    
    do j =1, jsize
       do i =1, isize
          !        print *,i,j
          if((mask(i,j).eq.999)) cycle
          !-----------------initialization---------------------------
          window=nint(oscale(i,j)*2.0*2.0)    !1.5 times radius
!          window=17
          mi=mod(window,2)
          if(mi.eq.0) window=window+1
          rad=(window-1)/2
          arad=(real(window-1))/2

          ii1=i-rad-1;jj1=j-rad-1
          if((ii1.lt.1).or.(jj1.lt.1)) cycle
          ii1=i+rad+1;jj1=j+rad+1                   
          if((ii1.gt.isize).or.(jj1.gt.jsize)) cycle 

          allocate(x1(window,window),y1(window,window))
          allocate(x2(window,window),y2(window,window))
          allocate(minx(window,window),miny(window,window))
          allocate(conv(window,window,4))

          !------------------ three x-directional scan and three y-directional scan per a window----------------
          bbp_max=-1
          do nt =1, n_templ
             do kk =1,2
                bbp(i,j)=0
                if(kk.eq.1) theta=(tmp_orient(nt)-orient(i,j))*dora
                if(kk.eq.2) theta=((tmp_orient(nt)-orient(i,j))-180)*dora

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
                enddo

                conv(:,:,:)=0
                !-------------------3% min- 97% max streching---------------------
                inum2=0              
                do jj=1,window
                   do ii=1,window                   
                      ai1=x2(ii,jj);aj1=y2(ii,jj)
                      if((ai1.lt.-0.51).or.(ai1.gt.0.51).or.(aj1.lt.-0.51).or.(aj1.gt.0.51)) cycle                    
                      inum2=1+inum2
                   enddo
                enddo

                allocate(list(inum2))

                do k=1,4
                   inum2=0 
                   do jj=1,window
                      do ii=1,window
                         ai1=x2(ii,jj);aj1=y2(ii,jj)
                         if((ai1.lt.-0.51).or.(ai1.gt.0.51).or.(aj1.lt.-0.51).or.(aj1.gt.0.51)) cycle  
                         inum2=1+inum2
                         !                    ii1=int(x1(ii,jj))+i;jj1=int(y1(ii,jj))+j
                         ii1=i+ii-rad-1;jj1=j+jj-rad-1                     
                         list(inum2)=data_i(ii1,jj1,k)
                      enddo
                   enddo

                   call heapsort2(inum2,list)
                   imax=nint(0.97*real(inum2))
                   imin=nint(0.03*real(inum2))
                   a=1/(list(imax)-list(imin))
                   b=list(imin)
                   do jj=1,window
                      do ii=1,window
                         ii1=i+ii-1-(window-1)/2;jj1=j+jj-1-(window-1)/2                  
                         conv(ii,jj,k)=a*(data_i(ii1,jj1,k)-b)
                         !a*(data_i(ii1,jj1,k)-b)
                      enddo
                   enddo
                   !                 print *,conv(:,2,k)
                enddo
                deallocate(list)
                !----------------------estimate BBPs-------------------------------
                inum=0
                minx(:,:)=999
                miny(:,:)=999
                minx2(:,:)=999
                miny2(:,:)=999                    
                dis(:,:)=10000.
                sum=0
                !        print *,i,j,step,window

                do jj=1,window
                   do ii=1,window
                      distance=1000000.
                      !                    ii1=int(x1(ii,jj))+i;jj1=int(y1(ii,jj))+j
                      ii1=i+ii-rad-1;jj1=j+jj-rad-1                    
                      ai1=x2(ii,jj);aj1=y2(ii,jj)
                      !                      print *,ii1,jj1,ii,jj,i,j

                      if((ii1.lt.1).or.(ii1.gt.isize).or.(jj1.lt.1).or.(jj1.gt.jsize)) cycle
                      if((ai1.lt.-0.51).or.(ai1.gt.0.51).or.(aj1.lt.-0.51).or.(aj1.gt.0.51)) cycle
                      sum=sum+1              
                      do jjj=1,jmax_tl,psize
                         do iii=1,imax_tl,psize
                         
                            ai2=x4(iii,jjj,nt);aj2=y4(iii,jjj,nt)
                            if((ai2.lt.-0.51).or.(ai2.gt.0.51).or.(aj2.lt.-0.51).or.(aj2.gt.0.51))  cycle

                            dis_xy=(dble(ai1-ai2))**2+(dble(aj1-aj2))**2
                            dis_rgb=0.
                            do k=1,band
                               dis_rgb=(dble(conv(ii,jj,k)-templ(iii,jjj,k,nt)))**2 + dis_rgb
                            enddo
                            adis=dis_rgb+dis_xy*disdep
                            if(adis.lt.dis(iii,jjj)) then
                               minx2(iii,jjj)=ii
                               miny2(iii,jjj)=jj
                               dis(iii,jjj)=adis
                            endif
                            if(adis.lt.distance) then
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

                      if((ii.eq.i3).and.(jj.eq.j3)) then
                         bbp(i,j)=bbp(i,j)+1
                      endif

                   enddo
                enddo
                
                sum_min=min(real(sum),real(imax_tl*jmax_tl/(psize**2)))
                bbp(i,j)=bbp(i,j)/sum_min-bbs_ref(int(sum_min))

                   
                if(bbp(i,j).ge.bbp_max) then
                   bbp_max=bbp(i,j)
                endif

             enddo
          enddo

          bbp(i,j)=bbp_max
          print *,bbp(i,j),i,j,oscale(i,j),orient(i,j),bbs_ref(int(sum_min))
          deallocate(x1,y1,x2,y2,minx,miny,conv)

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

  SUBROUTINE ARRANGE_TEMPLATE(tmpl2,isize2,jsize2,band2,n_templ,tmpl,isize,jsize,orient)
    IMPLICIT NONE
    INTEGER :: i,j,k,isize,jsize,band,ii,jj,i2,j2,l
    REAL(4):: xrate,yrate,xoffs,yoffs
  INTEGER(4),INTENT(IN) :: isize2,jsize2,band2,n_templ
  REAL(4),INTENT(IN) :: tmpl(isize,jsize,band2,n_templ)  
  REAL(4),INTENT(OUT) :: tmpl2(isize2,jsize2,band2,n_templ)
  character(512) :: infile
  REAL(4) :: gain(4),offset(4)
  REAL(4),ALLOCATABLE :: dummy(:,:,:)
  REAL(4) :: sum 
  REAL(4),INTENT(OUT) :: orient(n_templ)
  INTEGER(4),PARAMETER :: interval = 10
  REAL(4) :: hist(180/interval)
  REAL(4) :: dx,dy,pi,dora,theta,hismax,mu
  INTEGER :: itheta
  REAL(4) :: ndvi(isize2,jsize2,n_templ)  
!!$  
!!$  infile='template.tif'
!!$  infile=trim(dir)//infile  
!!$  call read_size_tif(infile,isize,jsize,band)
!!$  allocate(dummy(isize,jsize,band))
!!$  CALL open_tif(infile,dummy(:,:,:),isize,jsize,band,jsize)

  
  tmpl2(:,:,:,:)=0
  xrate=real(isize)/real(isize2)
  yrate=real(jsize)/real(jsize2)
  xoffs=(xrate+1.00000001)/2.
  yoffs=(yrate+1.00000001)/2.
  do l =1,n_templ
     do j =1,jsize2
        do i=1,isize2
           do k =1,band2           
           ii=nint(xoffs+real(i-1)*xrate)
           jj=nint(yoffs+real(j-1)*yrate)
           tmpl2(i,j,k,l)=tmpl(ii,jj,k,l)
        enddo
           ndvi(i,j,l)=(tmpl2(i,j,4,l)-tmpl2(i,j,3,l))/(tmpl2(i,j,4,l)+tmpl2(i,j,3,l))
        enddo
     enddo
  enddo
  
  pi=acos(-1.)
  dora=pi/180

  do l =1,n_templ
     hist(:)=0
     do jj=2,jsize2-1
        do ii=2,isize2-1
           i2=ii
           j2=jj
           dx=ndvi(ii+1,jj,l)-ndvi(ii-1,jj,l)
           dy=ndvi(ii,jj+1,l)-ndvi(ii,jj-1,l)
           mu=dx**2+dy**2           
           if(mu.eq.0) cycle
           theta=atan(dy/dx)/dora
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
           if(mu.le.90) orient(l)=180.-(mu+90)
           if(mu.ge.90) orient(l)=180.-(mu-90)
        endif
     enddo

  enddo

END SUBROUTINE ARRANGE_TEMPLATE

SUBROUTINE RESIZE(data_i,band2,jsize2,isize2,data_o,band,jsize,isize)
  IMPLICIT NONE
  INTEGER :: i,j,k,ii,jj
  REAL(4):: xrate,yrate,xoffs,yoffs
  INTEGER(4),INTENT(IN) :: isize2,jsize2,band2
  INTEGER(4),INTENT(IN) :: isize,jsize,band
  INTEGER(4), INTENT(OUT)  :: data_o(band,jsize,isize)
  REAL(4), INTENT(IN) :: data_i(band2,jsize2,isize2)  
  
  xrate=real(isize2)/real(isize)
  yrate=real(jsize2)/real(jsize)
  xoffs=(xrate+1.00000001)/2.
  yoffs=(yrate+1.00000001)/2.
     do j =1,jsize
        do i=1,isize
           do k =1,band           
           ii=nint(xoffs+real(i-1)*xrate)
           jj=nint(yoffs+real(j-1)*yrate)
           data_o(k,j,i)=data_i(k,jj,ii)
        enddo
     enddo
  enddo


END SUBROUTINE RESIZE

END MODULE Template_Matching_lr
