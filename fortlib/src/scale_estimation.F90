MODULE advanced_process
  USE Basic_processing
  USE similarity
  IMPLICIT NONE
  INTEGER, PARAMETER :: MAX_IO = 64

CONTAINS
  
 SUBROUTINE scale_selection(data_i,confidence_thres,scale_thres,data_o,band,jsize,isize)
  INTEGER :: i,j,k,l,i2,j2,kk
  INTEGER, PARAMETER :: band1=3
  INTEGER, PARAMETER :: band2=4
  INTEGER(4),INTENT(IN) :: isize,jsize,band
  REAL(4),INTENT(IN) :: confidence_thres,scale_thres
  REAL(4), INTENT(OUT)  :: data_o(2,jsize,isize)
  REAL(4), INTENT(IN) :: data_i(band,jsize,isize)
  REAL(4), ALLOCATABLE :: ndvi(:,:),data(:,:,:)
  INTEGER(4), ALLOCATABLE :: mask1(:,:),mask2(:,:),mask3(:,:)
  REAL(4), ALLOCATABLE :: octave(:,:,:),oscale(:,:)
  REAL(4), ALLOCATABLE :: diff(:,:,:),diff2(:,:,:)
  REAL(4), ALLOCATABLE :: dummy(:,:),confidence(:,:)
  INTEGER(4), ALLOCATABLE :: dummy2(:,:)
  INTEGER(1), ALLOCATABLE :: omax_scale(:,:,:)
  INTEGER(4) :: idirect,jdirect
  !------test data---
  REAL(4), ALLOCATABLE :: test(:,:,:)
  !---------
  INTEGER, PARAMETER :: octaves = 7
  INTEGER, PARAMETER :: interval = 2
  INTEGER, PARAMETER :: numpero = 5
  REAL, PARAMETER :: sigma=1
  REAL(4) :: sum,ksigma,k_int,diffmax,sf,diffmin,diff_dum,sum_conf,sum_conf2
  INTEGER :: mi,mj,wsize,dum_scale,disize,djsize,n_oc,interval2,first,sum_max
  INTEGER(4) :: isize2,jsize2,isize3,jsize3,conf_index,window,diffmax_scale

  !---------------INITIALIZATION-----------------------------
  allocate(ndvi(isize,jsize))  
  allocate(mask1(isize,jsize),mask3(isize,jsize))
  allocate(data(isize,jsize,band))

  CALL reshape_backforward_3D_REAL(data_i(1:band,1:jsize,1:isize),data(1:isize,1:jsize,1:band),band,jsize,isize)

  
  !--------------MAKE MASK DATA USING LAND COVER DATA-----------------------
  CALL Index(data(1:isize,1:jsize,1:band),ndvi(1:isize,1:jsize),isize,jsize,band,band1,band2,0)
  CALL Make_MASK_NODATA(ndvi(:,:),mask3(:,:),isize,jsize,999)
  mask1(:,:)=mask3(:,:)

!!$  do j =1,jsize
!!$     do i=1,isize
!!$        if((data(i,j,1).ge.0.99).and.(data(i,j,2).ge.0.99).and.(data(i,j,3).ge.0.99)) then
!!$           mask3(i,j)=999
!!$        endif
!!$     enddo
!!$  enddo

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
           endif
           if(diff2(i,j,k).lt.diffmin) then
              diffmin=diff2(i,j,k)
           endif           
        enddo
        if(diffmax.le.0.0) cycle
        oscale(i,j)=sigma*(k_int**real(dum_scale-1))*1.5
        confidence(i,j)=diffmax !(diffmax+diffmin)*(diffmax-diffmin)
     enddo
  enddo
  
  deallocate(diff,diff2,omax_scale)

  
  allocate(test(isize,jsize,2))
  test(:,:,1)=oscale(:,:)
  test(:,:,2)=confidence(:,:)
  CALL reshape_backforward_3D_REAL(test(:,:,:),data_o(:,:,:),isize,jsize,2)  

END SUBROUTINE scale_selection

  SUBROUTINE create_window_candidate(confidence0,oscale0,mask0,data_o,jsize,isize,confidence_thres,scale_thres)
  INTEGER :: i,j,k,l,ii,jj,kk,i2,j2,k2
  INTEGER(4),INTENT(IN) :: isize,jsize
  REAL(4),INTENT(IN) :: confidence_thres,scale_thres
  INTEGER(4), INTENT(OUT)  :: data_o(jsize,isize)
  REAL(4), INTENT(IN) :: confidence0(jsize,isize)
  INTEGER(4),INTENT(IN) :: mask0(jsize,isize)
  REAL(4), INTENT(IN) :: oscale0(jsize,isize)
  REAL(4) :: confidence(isize,jsize)
  INTEGER(4) :: mask(isize,jsize),mask2(isize,jsize)
  REAL(4) :: oscale(isize,jsize)  
  INTEGER(4) :: mi,window
  REAL(4) ::dist,overw
  !------test data---
  INTEGER(4), ALLOCATABLE :: test(:,:,:)
  INTEGER(4) :: ix,iy,sum
  
  CALL reshape_backforward_2D_INT(mask0,mask,jsize,isize)
  CALL reshape_backforward_2D_REAL(confidence0,confidence,jsize,isize)
  CALL reshape_backforward_2D_REAL(oscale0,oscale,jsize,isize)  
  !-------------select (mask) candidate windows------------------
  !-----1. confidence mask
  !-----2. mask unnecessary window using scale (overwrap) and confidence
  !-----3. scale mask (this mask must be put after the above mask 2)
  mask2(:,:)=mask(:,:)
  do j =1, jsize
     do i =1, isize
        if(mask(i,j).eq.999) cycle
        if(confidence(i,j).le.confidence_thres) then     
           mask(i,j)=999;mask2(i,j)=999
        endif
     enddo
  enddo
!!$!-----------------------confidence mask-------------------------------
!!$  do j =1, jsize
!!$     do i =1, isize
!!$        if(mask(i,j).eq.999) cycle
!!$
!!$        window=nint(oscale(i,j)*2)
!!$        mi=mod(window,2)
!!$        if(mi.eq.0) window=window+1
!!$        k=(window-1)/2
!!$
!!$        if((i+k.gt.isize).or.(i-k.lt.1)) then
!!$        mask(i,j)=999
!!$        cycle
!!$        endif
!!$        if((j+k.gt.jsize).or.(j-k.lt.1)) then
!!$        mask(i,j)=999
!!$        cycle
!!$        endif
!!$
!!$        do jj=1,window
!!$           do ii=1,window
!!$              i2=i+ii-1-k
!!$              j2=j+jj-1-k
!!$              if((i.eq.i2).and.(j.eq.j2)) cycle
!!$              
!!$              dist=sqrt(real((i-i2)**2+(j-j2)**2))
!!$              overw=real(oscale(i,j)/2.+oscale(i2,j2)/2.-dist)/real(oscale(i,j)/2.)
!!$
!!$              if(overw.lt.0.25) cycle
!!$              if(confidence(i,j).lt.confidence(i2,j2)) then
!!$                 mask(i,j)=999
!!$              elseif(confidence(i,j).gt.confidence(i2,j2)) then
!!$                 mask(i2,j2)=999
!!$              endif
!!$           enddo
!!$        enddo
!!$     enddo
!!$  enddo

!-----------------------confidence mask-------------------------------
  do j =1, jsize
     do i =1, isize
        if(mask(i,j).eq.999) cycle

        window=nint(oscale(i,j)*2)
        mi=mod(window,2)
        if(mi.eq.0) window=window+1
        k=(window-1)/2

        if((i+k.gt.isize).or.(i-k.lt.1)) then
        mask2(i,j)=999
        cycle
        endif
        if((j+k.gt.jsize).or.(j-k.lt.1)) then
        mask2(i,j)=999
        cycle
        endif

        do jj=1,window
           do ii=1,window
              i2=i+ii-1-k
              j2=j+jj-1-k
              if((i.eq.i2).and.(j.eq.j2)) cycle
              
              dist=sqrt(real((i-i2)**2+(j-j2)**2))
              overw=real(oscale(i,j)/2.+oscale(i2,j2)/2.-dist)/real(oscale(i,j)/2.)

              if(overw.lt.0.5) cycle
              if(confidence(i,j).lt.confidence(i2,j2)) then
                 mask2(i,j)=999;exit
                 exit
              endif
           enddo
          if(mask2(i,j).eq.999) exit
        enddo
     enddo
  enddo

!-----------------------confidence mask-------------------------------
  do j =1, jsize
     do i =1, isize
        if(mask2(i,j).eq.999) cycle

        window=nint(oscale(i,j)*2)
        mi=mod(window,2)
        if(mi.eq.0) window=window+1
        k=(window-1)/2

        do jj=1,window
           do ii=1,window
              i2=i+ii-1-k
              j2=j+jj-1-k
              if((i.eq.i2).and.(j.eq.j2)) cycle
              
              dist=sqrt(real((i-i2)**2+(j-j2)**2))
              overw=real(oscale(i,j)/2.+oscale(i2,j2)/2.-dist)/real(oscale(i,j)/2.)

              if(overw.lt.0.25) cycle
              if(confidence(i,j).lt.confidence(i2,j2)) then
                 mask2(i,j)=999
              elseif(confidence(i,j).gt.confidence(i2,j2)) then
                 mask2(i2,j2)=999
              endif
           enddo
        enddo
     enddo
  enddo
  
!------------------------------------------------------------
  do j =1, jsize
     do i =1, isize
        if(mask2(i,j).eq.999) cycle

        window=nint(oscale(i,j)*2)
        mi=mod(window,2)
        if(mi.eq.0) window=window+1
        k=(window-1)/2

        if((i+k.gt.isize).or.(i-k.lt.1)) cycle
        if((j+k.gt.jsize).or.(j-k.lt.1)) cycle

        ix=i
        iy=j
        sum=1
        mask2(i,j)=999
        
        do jj=1,window
           do ii=1,window
              i2=i+ii-1-k
              j2=j+jj-1-k
              if((i.eq.i2).and.(j.eq.j2)) cycle
              if(mask2(i2,j2).eq.999) cycle
              
             dist=sqrt(real((i-i2)**2+(j-j2)**2))
             overw=real(oscale(i,j)/2.+oscale(i2,j2)/2.-dist)/real(oscale(i,j)/2.)
              if(overw.lt.0.25) cycle
              mask2(i2,j2)=999
              ix=i2+ix
              iy=j2+iy
              sum=1+sum 
           enddo
        enddo

        mask2(int(real(ix)/real(sum)),int(real(iy)/real(sum)))=1
        
     enddo
  enddo
 
  do j =1, jsize
     do i =1, isize
        if(mask2(i,j).eq.999) cycle
         if(oscale(i,j).le.scale_thres) then
           mask2(i,j)=999
        endif
     enddo
  enddo

  CALL reshape_backforward_2D_INT(mask2(:,:),data_o(:,:),isize,jsize)  

END SUBROUTINE create_window_candidate


  SUBROUTINE orientation_selection(data_i,oscale0,mask0,data_o,band,jsize,isize)
    IMPLICIT NONE
  INTEGER :: i,j,k,ii,jj,kk,i2,j2,k2,i3,j3,k3,iii,jjj 
  INTEGER(4),INTENT(IN) :: isize, jsize, band
  REAL(4),INTENT(IN) :: data_i(band,jsize,isize)
  INTEGER(4),INTENT(IN) :: mask0(jsize,isize)
  REAL(4),INTENT(IN) :: oscale0(jsize,isize)
  REAL(4),INTENT(OUT):: data_o(jsize,isize)  
  INTEGER(4),PARAMETER :: interval = 10
  INTEGER(4),PARAMETER :: cal_n = 31
  INTEGER(4) :: mask(isize,jsize)
  REAL(4) :: oscale(isize,jsize)
  REAL(4),ALLOCATABLE :: data(:,:,:),ndvi(:,:),dum(:,:)
  real(4) :: tmp(0:cal_n+1,0:cal_n+1)
  REAL(4) :: hist(180/interval)
  REAL(4) :: dx,dy,pi,dora,theta,hismax,mu
  REAL(4) :: xrate,yrate,xoffs,yoffs
  INTEGER :: itheta,window,mi,step,window2
  REAL(4) :: sum
  REAL(4),ALLOCATABLE :: conv(:,:)
  INTEGER, PARAMETER :: band1=3
  INTEGER, PARAMETER :: band2=4
  
  
  pi=acos(-1.)
  dora=pi/180
  
  allocate(ndvi(isize,jsize),data(isize,jsize,band),dum(isize,jsize))
  CALL reshape_backforward_2D_INT(mask0,mask,jsize,isize)
  CALL reshape_backforward_2D_REAL(oscale0,oscale,jsize,isize)  
  CALL reshape_backforward_3D_REAL(data_i(1:band,1:jsize,1:isize),data(1:isize,1:jsize,1:band),band,jsize,isize)
  CALL Index(data(1:isize,1:jsize,1:band),ndvi(1:isize,1:jsize),isize,jsize,band,band1,band2,0)
  dum(:,:)=999.
  
  do j =1, jsize
     do i =1, isize
        if(mask(i,j).eq.999) cycle        
        window=nint(oscale(i,j)*2.0*1.5)    !1.5 times radius
        
        mi=mod(window,2)
        if(mi.eq.0) window=window+1
        k=(window-1)/2

        if((i+k+1.ge.isize).or.(i-k-1.le.1)) cycle
        if((j+k+1.ge.jsize).or.(j-k-1.le.1)) cycle        
        xrate=real(window)/real(cal_n)
        yrate=real(window)/real(cal_n)
        
        xoffs=real(i)-k-1
        yoffs=real(j)-k-1
        do j2=0,int(cal_n)+1
           do i2=0,int(cal_n)+1
           ii=nint(xoffs+real(i2-1)*xrate)
           jj=nint(yoffs+real(j2-1)*yrate)
           tmp(i2,j2)=ndvi(ii,jj)
        enddo
     enddo
     
     window=cal_n
        mi=mod(window,2)
        if(mi.eq.0) window=window+1
        k=(window-1)/2

        window2=5 !window/5+1
        mi=mod(window2,2)
        if(mi.eq.0) window2=window2+1
        k2=(window2-1)/2
        allocate(conv(0:window+1,0:window+1))
           do jj =0,window+1
              do ii =0,window+1
                 conv(ii,jj)=0
                 sum=0
                 do jjj = 1,window2
                    do iii =1,window2
                       i3=ii+iii-k2-1
                       j3=jj+jjj-k2-1                      
                       if((i3.gt.cal_n+1).or.(i3.lt.0)) cycle
                       if((j3.gt.cal_n+1).or.(j3.lt.0)) cycle
                       sum=1+sum
                       conv(ii,jj)=tmp(i3,j3)+conv(ii,jj)
                    enddo
                 enddo
                 conv(ii,jj)=conv(ii,jj)/sum
              enddo
           enddo

        hist(:)=0
        do jj=1,window
           do ii=1,window
              i2=ii-1-k
              j2=jj-1-k
              dx=conv(ii+1,jj)-conv(ii-1,jj)
              dy=conv(ii,jj+1)-conv(ii,jj-1)              
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
              if(mu.le.90) dum(i,j)=180.-(mu+90)
              if(mu.ge.90) dum(i,j)=180.-(mu-90)
           endif
        enddo
        deallocate(conv)
     enddo
  enddo

  CALL reshape_backforward_2D_REAL(dum(:,:),data_o(:,:),isize,jsize)  
  
END SUBROUTINE ORIENTATION_SELECTION

SUBROUTINE template_matching2(data_i,oscale0,confidence0,orient0,mask0,tmp_i,sim_thres1, &
     sim_thres2,sim_thres3,sim_thres4,scan,&
     data_o,isize,jsize,band,n_templ,band_t,jsize_t,isize_t)
  IMPLICIT NONE
  INTEGER(4),INTENT(IN) :: isize, jsize,band
  INTEGER(4),INTENT(IN) :: n_templ,band_t,jsize_t,isize_t
  REAL(4),INTENT(IN) :: sim_thres1,sim_thres2,sim_thres3,sim_thres4
  REAL(4),INTENT(OUT):: data_o(2,jsize,isize)
  REAL(4),INTENT(IN) :: data_i(band,jsize,isize)
  REAL(4),INTENT(IN) :: oscale0(jsize,isize)
  REAL(4),INTENT(IN) :: orient0(jsize,isize)
  INTEGER(4),INTENT(IN) :: mask0(jsize,isize)
  REAL(4),INTENT(IN) :: confidence0(jsize,isize)
  REAL(4), INTENT(IN) :: tmp_i(n_templ,band_t,jsize_t,isize_t)
  INTEGER(4),INTENT(IN) :: scan !must be a odd number

  REAL(4) :: data(isize,jsize,band)   
  REAL(4) :: tmp0(isize_t,jsize_t,band,n_templ)
  REAL(4) :: oscale(isize,jsize)
  REAL(4) :: orient(isize,jsize)
  INTEGER(4) :: mask(isize,jsize)
  REAL(4) :: confidence(isize,jsize)

  REAL(4) :: tmp_orient(n_templ)
  INTEGER(4), PARAMETER :: imax_tl=17 !17
  INTEGER(4), PARAMETER :: jmax_tl=17 !17
  REAL(4) :: templ(imax_tl,jmax_tl,band,n_templ)
  REAL(4) :: bright_t(imax_tl,jmax_tl,1,n_templ)
  REAL(4) :: ndvi_t(imax_tl,jmax_tl,1,n_templ)
  REAL(4) :: gsi_t(imax_tl,jmax_tl,1,n_templ)
  
  INTEGER :: i,j,ii,jj,k,kk,iii,jjj,i2,j2,i3,j3,i4,j4,scan_i,scan_j
  INTEGER(4) :: ii1,ii2,jj1,jj2,window
  INTEGER(4) :: imin,imax,jmin,jmax,lwindow
  INTEGER(4) :: inum,lwindow_tl
  INTEGER(4) :: mi,mj,isize2,jsize2,rad,nt,step

  REAL(4) :: sim(isize,jsize)
  REAL(4) :: sim2(isize,jsize)
  REAL(4) :: sim3(isize,jsize)
  REAL(4) :: sim4(isize,jsize)    
  REAL(4) :: sim_scan(scan,scan)
  REAL(4) :: sim2_scan(scan,scan)
  REAL(4) :: sim3_scan(scan,scan)
  REAL(4) :: sim4_scan(scan,scan) 
  REAL(4) :: sim0,sim02,sim03,sim04
  REAL(4) :: sim_thres22,sim_thres32,sim_thres42

  REAL(4), ALLOCATABLE :: test(:,:,:)
  REAL(4),ALLOCATABLE :: forshape(:,:) 
  REAL(4),ALLOCATABLE :: x1(:,:),x2(:,:),x3(:,:,:),x4(:,:,:)
  REAL(4),ALLOCATABLE :: y1(:,:),y2(:,:),y3(:,:,:),y4(:,:,:)
  REAL(4),ALLOCATABLE :: ndvi(:,:,:),conv0(:,:,:),conv(:,:,:)
  REAL(4),ALLOCATABLE :: brightness(:,:,:),gsi(:,:,:)
  
  REAL(4) :: ai1,ai2,aj1,aj2,pi,dora,theta
  REAL(4) :: sim_max,sim_dum
  REAL(4) :: sim2_argmax,sim3_argmax,sim4_argmax
  real(4) :: sum,arad
  real(4) :: ndvi2,si2
  real(4) :: scan2
  INTEGER(4) :: conf_index,mai
  CHARACTER(512):: file
  CALL reshape_backforward_2D_INT(mask0,mask,jsize,isize)
  CALL reshape_backforward_2D_REAL(oscale0,oscale,jsize,isize)
  CALL reshape_backforward_2D_REAL(orient0,orient,jsize,isize)
  CALL reshape_backforward_2D_REAL(confidence0,confidence,jsize,isize)   
  CALL reshape_backforward_4D_REAL(tmp_i,tmp0,n_templ,band_t,jsize_t,isize_t)
  CALL ARRANGE_TEMPLATE(templ,imax_tl,jmax_tl,band_t,n_templ,tmp0,isize_t,jsize_t,tmp_orient)
  CALL reshape_backforward_3D_REAL(data_i(1:band,1:jsize,1:isize),data(1:isize,1:jsize,1:band),band,jsize,isize)  
  print *,"template orientation =", tmp_orient
  pi=acos(-1.)
  dora=pi/180.
  sim(:,:)=0
  sim2(:,:)=0

  if(sim_thres2.ge.0) then
  sim_thres22=sim_thres2
  else
  sim_thres22=0.35 * sim_thres2
  endif

  if(sim_thres4.ge.0) then
  sim_thres42=sim_thres4
  else
  sim_thres42=0.5 * sim_thres4
  endif

  if(sim_thres3.ge.0) then
  sim_thres32= 0.95 * sim_thres3 + 0.05
  else
  sim_thres32= 1.05 * sim_thres3 + 0.05
  endif

  !-----------calc distance (-0.5 - 0.5). must set theta=0 (;template image projection base)------
  allocate(x3(imax_tl,jmax_tl,n_templ),x4(imax_tl,jmax_tl,n_templ))
  allocate(y3(imax_tl,jmax_tl,n_templ),y4(imax_tl,jmax_tl,n_templ))
  theta=0
  do nt=1,n_templ
     call calc_rot_loc(x3(:,:,nt),y3(:,:,nt),x4(:,:,nt),y4(:,:,nt),imax_tl,jmax_tl,theta)
  enddo

  !-------------------normalization or standardization---------------------   
  do nt=1,n_templ
     ndvi_t(:,:,1,nt) = (templ(:,:,4,nt)-templ(:,:,3,nt))/(templ(:,:,4,nt)+templ(:,:,3,nt))
     gsi_t(:,:,1,nt) = (templ(:,:,3,nt)-templ(:,:,1,nt))/(templ(:,:,3,nt)+templ(:,:,2,nt)+templ(:,:,1,nt))
     bright_t(:,:,1,nt) = (templ(:,:,1,nt)+templ(:,:,2,nt)+templ(:,:,3,nt))/3.
     call normalization_local(templ,imax_tl,jmax_tl,band)
!     call standardization_local(templ,imax_tl,jmax_tl,band)
  enddo

  !-------------------------------------------------------------------------------
  !----------------------------- main loop ---------------------------------------
  !-------------------------------------------------------------------------------   
  scan2=(scan-1)/2
  lwindow_tl=imax_tl + int((imax_tl-1)/2./4.)*scan2*2
  allocate(x1(imax_tl,jmax_tl),y1(imax_tl,jmax_tl))
  allocate(x2(imax_tl,jmax_tl),y2(imax_tl,jmax_tl))
  allocate(conv0(lwindow_tl,lwindow_tl,4))
  allocate(conv(imax_tl,jmax_tl,4))
  allocate(brightness(imax_tl,jmax_tl,1))
  allocate(ndvi(imax_tl,jmax_tl,1))
  allocate(gsi(imax_tl,jmax_tl,1))  
  
  do j =1,jsize!474,474!1, jsize
     print *,j,jsize
     do i =1,isize!4035, 4035 !1, isize
        
        if((mask(i,j).eq.999.).or.(orient(i,j).eq.999)) cycle
        !-----------------initialization---------------------------
        window=nint(oscale(i,j)*2.0*1.5)   !4.0 times windows 
        mi=mod(window,2)
        if(mi.eq.0) window=window+1
        arad=(real(window-1))/2
        rad=(window-1)/2

        imin=i+int(arad/4.)*(-scan2)-rad
        imax=i+int(arad/4.)*(scan2)+rad
        jmin=j+int(arad/4.)*(-scan2)-rad   
        jmax=j+int(arad/4.)*(scan2)+rad
        lwindow=imax-imin+1
        
!        print *,oscale(i,j),confidence(i,j),orient(i,j)
!        print *,imin,imax,jmin,jmax,lwindow,lwindow_tl,window
        if((imin.lt.1).or.(imax.gt.isize).or.(jmin.lt.1).or.(jmax.gt.jsize)) cycle        !must modify
        
        mai=0
        do jj = jmin,jmax
        do ii = imin,imax
            if((data(ii,jj,1).le.0.01).and.(data(ii,jj,2).le.0.01).and.(data(ii,jj,3).le.0.01)) then
            mai=1+mai
            endif
        enddo
        enddo
        
        if(real(mai)/real(window*window).ge.0.2) cycle
        
        call gaussian_resize(data(imin:imax,jmin:jmax,:),conv0, &
             lwindow,lwindow,lwindow_tl,lwindow_tl,band)

        !------------------ three x- and y-directional scan per a window----------------
        do scan_j=1,scan
           do scan_i=1,scan
              imin=1 + int((imax_tl-1)/2./4.)*(scan_i-1)
              imax=imax_tl + int((imax_tl-1)/2./4.)*(scan_i-1)
              jmin=1 + int((jmax_tl-1)/2./4.)*(scan_j-1)
              jmax=jmax_tl + int((jmax_tl-1)/2./4.)*(scan_j-1)

              i4=i+int(arad/4.)*(scan_i-scan2-1)
              j4=j+int(arad/4.)*(scan_j-scan2-1)

              conv=conv0(imin:imax,jmin:jmax,:)
              brightness(:,:,1)=(conv(:,:,1)+conv(:,:,2)+conv(:,:,3))/3.
              ndvi(:,:,1)=(conv(:,:,4)-conv(:,:,3))/(conv(:,:,4)+conv(:,:,3))
              gsi(:,:,1)=(conv(:,:,3)-conv(:,:,1))/(conv(:,:,3)+conv(:,:,2)+conv(:,:,1))  
              
              call normalization_local(conv,imax_tl,jmax_tl,band)
              !-----------------------SIMs calculation for a window-------------
              sim_max=0
              sim2_argmax=0
              do nt =1, n_templ

                 do kk =1,2               !kk = a  loop for inverse rotations
                    sim(i4,j4)=0                   
                    if(kk.eq.1) theta=(tmp_orient(nt)-orient(i,j))*dora
                    if(kk.eq.2) theta=((tmp_orient(nt)-orient(i,j))-180)*dora

                    call calc_rot_loc(x1,y1,x2,y2,imax_tl,jmax_tl,theta)

                    !----------------------estimate SIMs-------------------------------
                    call calculate_bbs(conv,x2,y2,templ(:,:,:,nt),x4,y4,imax_tl,jmax_tl,band,1,sim0)
                    call calculate_zncc(brightness(:,:,:),bright_t(:,:,:,nt),imax_tl,jmax_tl,1,sim02)  
                    call calculate_zncc(ndvi(:,:,:),ndvi_t(:,:,:,nt),imax_tl,jmax_tl,1,sim03)
                    call calculate_zncc(gsi(:,:,:),gsi_t(:,:,:,nt),imax_tl,jmax_tl,1,sim04)                    
!                    call calculate_zncc(conv,templ(:,:,:,nt),imax_tl,jmax_tl,band,sim0)
                    if(sim0.ge.sim_max) then
                       sim_max=sim0
                       sim2_argmax = sim02
                       sim3_argmax = sim03
                       sim4_argmax = sim04                       
                    endif

                 enddo
              enddo

              sim_scan(scan_i,scan_j)=sim_max
              sim2_scan(scan_i,scan_j)=sim2_argmax
              sim3_scan(scan_i,scan_j)=sim3_argmax
              sim4_scan(scan_i,scan_j)=sim4_argmax             
           enddo
        enddo

         sim_max=-999
         do scan_j =1,scan
            do scan_i=1,scan
               i2=i+int(arad/4.)*(scan_i-scan2-1)
               j2=j+int(arad/4.)*(scan_j-scan2-1)
               if(sim_scan(scan_i,scan_j).gt.sim_max) then
                  i3=i2
                  j3=j2
                  sim_max=sim_scan(scan_i,scan_j)
                  sim2_argmax=sim2_scan(scan_i,scan_j)
                  sim3_argmax=sim3_scan(scan_i,scan_j)
                  sim4_argmax=sim4_scan(scan_i,scan_j)                   
               endif
            enddo
         enddo
         mask(i,j)=999
         oscale(i3,j3)=oscale(i,j)
         mask(i3,j3)=1
         sim(i3,j3)=sim_max
         sim2(i3,j3)=sim2_argmax
         sim3(i3,j3)=sim3_argmax
         sim4(i3,j3)=sim4_argmax           
!        print *,sim_scan(:,:)
!        print *,sim(i3,j3),i,j,oscale(i3,j3),window,orient(i3,j3)

     enddo
  enddo
  
  deallocate(x1,y1,x2,y2,conv,brightness,ndvi,gsi)

  allocate(forshape(isize,jsize))
  forshape(:,:)=0
  
  do j =1,jsize
     do i =1,isize
        if(mask(i,j).eq.999) cycle
        if(sim(i,j).lt.sim_thres1) cycle        
        if(sim2(i,j).lt.sim_thres22) cycle
        if(sim3(i,j).lt.sim_thres32) cycle
        if(sim4(i,j).lt.sim_thres42) cycle
        if((oscale(i,j).ge.60).and.(sim2(i,j).le.-0.2)) cycle        
!        if(log_reg_classifier(sim(i,j),sim2(i,j),sim3(i,j),sim4(i,j)).lt.0.5) cycle
        window=int(oscale(i,j)*1.5)
        mi=mod(window,2)
        if(mi.eq.0) window=window+1
        k=(window-1)/2
        conf_index=1
        do jj =1,window
           do ii =1,window
              i2=i+ii-k-1-2
              j2=j+jj-k-1-2
              if((i2.lt.1).or.(i2.gt.isize).or.(j2.lt.1).or.(j2.gt.jsize)) cycle
              if(sim(i2,j2).lt.sim_thres1) cycle
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

  allocate(test(isize,jsize,2))

!!$  test(:,:,1)=sim(:,:)
!!$  test(:,:,2)=forshape(:,:)

  test(:,:,1)=sim4(:,:)
  test(:,:,2)=forshape(:,:)
!  test(:,:,3)=sim2(:,:)
!  test(:,:,4)=sim3(:,:)  
!  print *,"hi"
  CALL reshape_backforward_3D_REAL(test(:,:,:),data_o(:,:,:),isize,jsize,2)  

  
 END SUBROUTINE template_matching2
 

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

!  tmpl2(:,:,:,:)=0
!  print *,"start filter"
!  do l =1,n_templ
!  call gaussian_resize(tmpl(:,:,:,l),tmpl2(:,:,:,l),isize,jsize,isize2,jsize2,band2)
!  enddo
!  print *,"finish filter"
  tmpl2(:,:,:,:) = tmpl

  xrate=real(isize)/real(isize2)
  yrate=real(jsize)/real(jsize2)
  xoffs=(xrate+1)/2.
  yoffs=(yrate+1)/2.
  do l =1,n_templ
     do j =1,jsize2
        do i=1,isize2
!!$           do k =1,band2           
!!$           ii=nint(xoffs+real(i-1)*xrate)
!!$           jj=nint(yoffs+real(j-1)*yrate)
!!$           tmpl2(i,j,k,l)=tmpl(ii,jj,k,l)
!!$        enddo
           ndvi(i,j,l)=(tmpl2(i,j,4,l)-tmpl2(i,j,3,l))/(tmpl2(i,j,4,l)+tmpl2(i,j,3,l))
        enddo
     enddo
  enddo

!!$  do j = 1, jsize2
!!$     print *,ndvi(:,j,1)
!!$  enddo
!!$  print *,n_templ
  
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
!     print *, "hist"
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
  REAL(4), INTENT(OUT)  :: data_o(band,jsize,isize)
  REAL(4), INTENT(IN) :: data_i(band2,jsize2,isize2)  
  REAL(4) :: data_i2(isize2,jsize2,band2)
  REAL(4) :: data_o2(isize,jsize,band)
  
!  xrate=real(isize2)/real(isize)
!  yrate=real(jsize2)/real(jsize)
!  xoffs=(xrate+1)/2.
!  yoffs=(yrate+1)/2.
!     do j =1,jsize
!        do i=1,isize
!           do k =1,band           
!           ii=nint(xoffs+real(i-1)*xrate)
!           jj=nint(yoffs+real(j-1)*yrate)
!           data_o(k,j,i)=data_i(k,jj,ii)
!        enddo
!     enddo
!  enddo
  CALL reshape_backforward_3D_REAL(data_i,data_i2,band2,jsize2,isize2)  
  call gaussian_resize(data_i2(:,:,:),data_o2(:,:,:),isize2,jsize2,isize,jsize,band2)
  CALL reshape_backforward_3D_REAL(data_o2,data_o,isize,jsize,band)  

END SUBROUTINE RESIZE

  SUBROUTINE calc_histgram_parameter(data_i,bin_width,bin_number,hist_par,hmax,hmode,pdf,band,jsize,isize)
  INTEGER(4) :: i, j,hi,inum,hsize
  INTEGER, PARAMETER :: band1=3
  INTEGER, PARAMETER :: band2=4
  INTEGER(4),INTENT(IN) :: isize,jsize,band,bin_number
  REAL(4),INTENT(IN) :: bin_width,hist_par
  REAL(4), INTENT(IN) :: data_i(band,jsize,isize)
  REAL(4), INTENT(OUT) :: hmax,hmode
  REAL(4) :: alist(isize*jsize)
  REAL(4), ALLOCATABLE :: list(:)
  INTEGER(4) :: smask(isize,jsize),smask2(isize,jsize)
  INTEGER(4), ALLOCATABLE :: morp1(:,:),morp2(:,:),bin(:,:)
  REAL(4), ALLOCATABLE :: ndvi(:,:),data(:,:,:)
  INTEGER(4), ALLOCATABLE :: mask1(:,:),mask2(:,:)

  REAL(4) :: amax
  REAL(4), ALLOCATABLE :: histgram(:)
  REAL(4), INTENT(OUT) :: pdf(2, bin_number)
  
  !---------------INITIALIZATION-----------------------------
  allocate(ndvi(isize,jsize),bin(isize,jsize))
  allocate(mask1(isize,jsize),mask2(isize,jsize))
  allocate(data(isize,jsize,band))

  CALL reshape_backforward_3D_REAL(data_i(1:band,1:jsize,1:isize),data(1:isize,1:jsize,1:band),band,jsize,isize)
  CALL Make_Mask_FORLANDSLIDE(data(:,:,4:5),smask(:,:),isize,jsize)
  
  !--------------MAKE MASK DATA USING LAND COVER DATA-----------------------
  CALL Index(data(1:isize,1:jsize,1:band),ndvi(1:isize,1:jsize),isize,jsize,band,band1,band2,0)
     
  inum=0
  amax=0
   do j = 1, jsize
      do i = 1, isize
         if((smask(i,j).eq.1).and.(ndvi(i,j).ne.999)) then
         inum=1+inum
         alist(inum)=ndvi(i,j)
         if(ndvi(i,j).ge.amax) amax=ndvi(i,j)
      endif
      enddo
   enddo
!   print *, amax
   allocate(list(inum))
   list(1:inum)=alist(1:inum)
   call heapsort2(inum,list)
   hmax=list(nint(hist_par*real(inum)))
   
  hsize=nint(1/bin_width)
  allocate(histgram(1:hsize))
  call calc_histgram(list(1:inum),bin_width,histgram,inum,hsize)
  
  amax=0
  do i = 1, hsize
     pdf(1,i)=real(bin_width*i)+real(bin_width/2.)
     pdf(2,i)=histgram(i)
     if(histgram(i).ge.amax) then
        amax=histgram(i)
        hi=i      
     endif
  enddo

  hmode=bin_width*real(hi) !-bin_width/2.

END SUBROUTINE calc_histgram_parameter



SUBROUTINE calc_histgram(list,bin,histgram,isize,hsize)
  INTEGER,INTENT(IN) :: isize,hsize
  REAL(4),INTENT(IN) :: bin
  REAL(4),INTENT(IN) :: list(isize)
  REAL(4),INTENT(OUT) :: histgram(hsize)
  REAL(4) :: dum(0:hsize)
  INTEGER :: i,j
  INTEGER :: inum

 dum(:)=0
 inum=0
 do i = 1, isize
    if(list(i).ge.bin*inum+bin/2.) then
       dum(inum)=real(i)
!       print *,dum(inum),inum,list(i)      
       inum=inum+1
    endif
 enddo

 do i =1,inum
    histgram(i)=dum(i)-dum(i-1)
    if(histgram(i).le.0) histgram(i)=0.
!    print *,histgram(i)
 enddo
 
ENDSUBROUTINE calc_histgram
END MODULE advanced_process
