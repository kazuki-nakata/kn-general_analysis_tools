MODULE similarity
  !USE Basic_processing
  IMPLICIT NONE

CONTAINS
  SUBROUTINE calculate_ssd(image,image2,isize,jsize,band,nst,similarity)
    IMPLICIT NONE    
    INTEGER(4) :: i, j, k
    INTEGER(4), INTENT(IN) :: isize,jsize,band,nst
    REAL(4), INTENT(OUT)  :: similarity
    REAL(4), INTENT(IN) :: image(isize,jsize,band)
    REAL(4), INTENT(IN) :: image2(isize,jsize,band)

    similarity = 0
    do j =1, jsize, nst
       do i =1, isize, nst
          do k =1, band
             similarity = similarity + (image(i,j,k) - image2(i,j,k))**2
          enddo
       enddo
    enddo

  END SUBROUTINE calculate_ssd

  SUBROUTINE calculate_sad(image,image2,isize,jsize,band, nst,similarity)
    IMPLICIT NONE    
    INTEGER(4) :: i, j, k
    INTEGER(4), INTENT(IN) :: isize,jsize,band, nst   
    REAL(4), INTENT(OUT)  :: similarity
    REAL(4), INTENT(IN) :: image(isize,jsize,band)
    REAL(4), INTENT(IN) :: image2(isize,jsize,band)

    similarity = 0
    do j =1, jsize, nst
       do i =1, isize, nst
          do k =1, band
             similarity = similarity + ABS(image(i,j,k) - image2(i,j,k))
          enddo
       enddo
    enddo

  END SUBROUTINE calculate_sad

  SUBROUTINE calculate_ncc(image,image2,isize,jsize,band, nst,similarity)
    IMPLICIT NONE    
    INTEGER(4) :: i, j, k
    INTEGER(4), INTENT(IN) :: isize,jsize,band, nst
    REAL(4), INTENT(OUT)  :: similarity
    REAL(4), INTENT(IN) :: image(isize,jsize,band)
    REAL(4), INTENT(IN) :: image2(isize,jsize,band)
    REAL(4) :: numerator, denominator, denominator1, denominator2

    numerator = 0
    denominator1 = 0
    denominator2 = 0
    do j =1, jsize, nst
       do i =1, isize, nst
          do k =1, band
             numerator = numerator + image(i,j,k)*image2(i,j,k)
             denominator1 = denominator1 + image(i,j,k)**2
             denominator2 = denominator2 + image2(i,j,k)**2
          enddo
       enddo
    enddo

    denominator = sqrt(denominator1)*sqrt(denominator2)

    if (denominator == 0) then
       similarity = 0
    else
       similarity = numerator / denominator
    endif

  END SUBROUTINE calculate_ncc

  SUBROUTINE calculate_zncc(image,image2,isize,jsize,band,nst,similarity)
    IMPLICIT NONE    
    INTEGER(4) :: i, j, k
    REAL(8) :: x1,y1,x2,y2,xy
    REAL(4), INTENT(OUT)  :: similarity
    INTEGER(4), INTENT(IN) :: isize,jsize,band, nst
    REAL(4), INTENT(IN) :: image(isize,jsize,band)
    REAL(4), INTENT(IN) :: image2(isize,jsize,band)
   !  REAL(4) :: diff(isize,jsize,band)
   !  REAL(4) :: diff2(isize,jsize,band)
    REAL(4) :: numerator, denominator, denominator1, denominator2,sum
    REAL(4) :: ave_image1(band), ave_image2(band) 

   sum=real(isize*jsize)
   !  ave_image1(:)=0.
   !  ave_image2(:)=0.
   !  do j = 1,jsize
   !  do i = 1, isize
   !  do k =1, band
   !   ave_image1(k)=image(i,j,k)/sum+ave_image1(k)
   !   ave_image2(k)=image2(i,j,k)/sum+ave_image2(k)
   !  enddo
   !  enddo
   !  enddo

   !  do k =1,band
   !  do j =1,jsize
   !  do i =1,isize
   !  diff(i,j,k) = image(i,j,k) - ave_image1(k)
   !  diff2(i,j,k) = image2(i,j,k) - ave_image2(k)    
   !  enddo
   !  enddo
   !  enddo
    

   !  numerator = 0
   !  denominator1 = 0
   !  denominator2 = 0
   !  do j =1, jsize
   !     do i =1, isize
   !        do k =1, band
   !           numerator = numerator + diff(i,j,k)*diff2(i,j,k)
   !           denominator1 = denominator1 + diff(i,j,k)**2
   !           denominator2 = denominator2 + diff2(i,j,k)**2
   !        enddo
   !     enddo
   !  enddo
    sum=0
    x1=0; x2=0; y1=0; y2=0; xy=0
    do j =1, jsize, nst
       do i =1, isize, nst
          do k =1, band
              sum=1+sum
              x1=image(i,j,k)+x1
              x2=image(i,j,k)**2+x2
              y1=image2(i,j,k)+y1
              y2=image2(i,j,k)**2+y2
              xy=image(i,j,k)*image2(i,j,k)+xy
          enddo
       enddo
    enddo

    denominator1=real(x2)/sum-(real(x1)/sum)**2
    denominator2=real(y2)/sum-(real(y1)/sum)**2
    numerator=real(xy)/sum-real(x1*y1)/(sum**2)

    denominator = sqrt(denominator1*denominator2)

    if (denominator == 0) then
       similarity = 0
    else
       similarity = numerator / denominator
    endif

  END SUBROUTINE calculate_zncc

  SUBROUTINE calculate_bbs(image,x,y,image2,x2,y2,isize,jsize,band,nst,disdep,sim)
    IMPLICIT NONE    
    INTEGER(4) :: i,j,k,i2,j2,i3,j3,i4,j4,ii,jj
    INTEGER(4), INTENT(IN) :: isize,jsize,band,nst
    REAL(4), INTENT(OUT)  :: sim
    REAL(4), INTENT(IN) :: x(isize,jsize),y(isize,jsize)
    REAL(4), INTENT(IN) :: x2(isize,jsize),y2(isize,jsize)
    REAL(4), INTENT(IN) :: image(isize,jsize,band)
    REAL(4), INTENT(IN) :: image2(isize,jsize,band)
    REAL(4), INTENT(IN) :: disdep
    REAL(4) :: ai1,aj1,ai2,aj2
    REAL(8) :: dis_xy,dis_rgb,adis,distance
    INTEGER(4) :: minx(isize,jsize),miny(isize,jsize)
    INTEGER(4) :: minx2(isize,jsize),miny2(isize,jsize)
    REAL(8) :: dis(isize,jsize)
    INTEGER(4) :: sum1,sum2
   !  REAL(8),PARAMETER :: disdep=0.5 !dependency factor of pixels distance
    
    sum1 = 0
    sum2 = isize*jsize/nst/nst
    minx(:,:)=999
    miny(:,:)=999
    minx2(:,:)=999
    miny2(:,:)=999                    
    dis(:,:)=10000.    
    
    
    do j=1,jsize, nst
       do i=1,isize, nst
          distance=1000000.
          ai1=x(i,j);aj1=y(i,j)
          if((ai1.lt.-0.51).or.(ai1.gt.0.51).or.(aj1.lt.-0.51).or.(aj1.gt.0.51)) cycle
          sum1=sum1+1              
          do jj=1,jsize, nst
             do ii=1,isize, nst
                ai2=x2(ii,jj);aj2=y2(ii,jj)
                if((ai2.lt.-0.51).or.(ai2.gt.0.51).or.(aj2.lt.-0.51).or.(aj2.gt.0.51)) cycle
                dis_xy=dble(ai1-ai2)**2+dble(aj1-aj2)**2
                dis_rgb=0
                do k=1,band
                   dis_rgb=dble(image(i,j,k)-image2(ii,jj,k))**2 + dis_rgb
                enddo
!                print *,dis_rgb,dis_xy
                adis=dis_rgb+dis_xy*disdep
                if(adis.le.dis(ii,jj)) then
                   minx2(ii,jj)=i
                   miny2(ii,jj)=j
                   dis(ii,jj)=adis
                endif
                if(adis.le.distance) then
                   minx(i,j) = ii
                   miny(i,j) = jj
                   distance=adis
                endif
             enddo
          enddo
       enddo
    enddo

!    print *,minx2(:,:)
!    print *,miny2(:,:)

    sim = 0
    do j=1,jsize
       do i=1,isize
!          print *,i,j
          if(minx(i,j).eq.999) cycle
          i2=minx(i,j)
          j2=miny(i,j)
!          print *,i2,j2,'2'
          i3=minx2(i2,j2)
          j3=miny2(i2,j2)
!          print *,i3,j3,'3'
          if((i.eq.i3).and.(j.eq.j3)) then
             sim=sim+1
!             print *,sim
           endif
       enddo
    enddo

    sim=sim/min(real(sum1),real(sum2))
!    print *, sim,"hi2"
 END SUBROUTINE calculate_bbs

SUBROUTINE calc_rot_loc(x2,y2,imax,jmax,theta)
  INTEGER(4),INTENT(IN) :: imax,jmax
  REAL(4),INTENT(IN) :: theta
  REAL(4) :: x(imax,jmax)
  REAL(4) :: y(imax,jmax)
  REAL(4),INTENT(OUT) :: x2(imax,jmax)
  REAL(4),INTENT(OUT) :: y2(imax,jmax)
  INTEGER(4) :: kx,ky,i,j,i2,j2
  REAL(4) :: maxdis,ai2,aj2

  kx=(imax-1)/2
  ky=(jmax-1)/2
  maxdis=sqrt(real(kx)**2+real(ky)**2)  
  do j = 1, jmax
     do i = 1, imax
        !           i2=i-kx-1
        !           j2=j-ky-1
        ai2=real(i)-real(kx)-1
        aj2=real(j)-real(ky)-1
        x(i,j)=ai2*cos(theta)+aj2*sin(theta)
        y(i,j)=aj2*cos(theta)-ai2*sin(theta)
        x2(i,j)=(ai2*cos(theta)+aj2*sin(theta))/real(imax-1)
        y2(i,j)=(aj2*cos(theta)-ai2*sin(theta))/real(jmax-1)
     enddo
  enddo
ENDSUBROUTINE calc_rot_loc

! SUBROUTINE normalization_local(data,imax,jmax,band)
!   INTEGER(4),INTENT(IN) :: imax,jmax,band
!   REAL(4),INTENT(INOUT) :: data(imax,jmax,band)
!   INTEGER(4) :: i,j,k,ii,jj,inum,imax2,imin2
!   REAL(4) :: a,b
!   REAL(4) :: list(imax*jmax)

! !!$  temp_sum=0     
! !!$    do jj = 1, jmax_tl
! !!$       do ii = 1, imax_tl
! !!$          if((x4(ii,jj,1).lt.-0.51).or.(x4(ii,jj,1).gt.0.51).or. &
! !!$               (y4(ii,jj,1).lt.-0.51).or.(y4(ii,jj,1).gt.0.51)) cycle
! !!$          temp_sum=temp_sum+1
! !!$       enddo
! !!$    enddo

!   do k =1,band
!      inum=0         
!      do j = 1, jmax
!         do i = 1, imax
!            inum=1+inum
!            list(inum)=data(i,j,k)
!         enddo
!      enddo

!      call heapsort2(inum,list)
!      imax2=nint(0.97*real(inum))
!      imin2=nint(0.03*real(inum))

!      a=1/(list(imax2)-list(imin2))

!      b=list(imin2)

!      data(:,:,k)=a*(data(:,:,k)-b)

!   enddo

! ENDSUBROUTINE normalization_local

SUBROUTINE standardization_local(data,imax,jmax,band)
  INTEGER(4),INTENT(IN) :: imax,jmax,band
  REAL(4),INTENT(INOUT) :: data(imax,jmax,band)
  INTEGER(4) :: i,j,k,ii,jj,inum,imax2,imin2
  REAL(4) :: sum,a,b
  REAL(4) :: ave(band),sigma(band)

  ave(:)=0
  sum=real(imax*jmax)
  do k =1,band
     do j = 1, jmax
        do i = 1, imax
           ave(k)=data(i,j,k)/sum+ave(k)
        enddo
     enddo
  enddo

  sigma(:)=0
  do k =1,band
     do j = 1, jmax
        do i = 1, imax
           sigma(k)=sigma(k) + (data(i,j,k)-ave(k))**2
        enddo
     enddo
  enddo

  sigma(:) = sqrt(sigma(:)/real(imax*jmax-1))

  do k =1,band
     data(:,:,k)=(data(:,:,k)-ave(k))/sigma(k)
  enddo

ENDSUBROUTINE standardization_local

SUBROUTINE gaussian_resize(input,output,isize,jsize,isize2,jsize2,band)
  IMPLICIT NONE
  INTEGER :: i,j,k,i2,j2,ii,jj
  INTEGER(4),INTENT(IN) :: isize,jsize,isize2,jsize2,band
  REAL(4),INTENT(OUT):: output(isize2,jsize2,band)
  REAL(4),INTENT(IN) :: input(isize,jsize,band)
  REAL(4) :: kernel
  INTEGER(4) :: wsize,reduce
  REAL(4) :: sigma,sum,sf,alpha
  REAL(4) :: mx,my,dx,dy,arad,xrate,yrate
  REAL(4) :: pi,sigma2


  output(:,:,:)=0
  xrate=real(isize)/real(isize2)
  yrate=real(jsize)/real(jsize2)

  wsize = int(xrate)+1
  sigma=(xrate+1)/6.
  arad = real(wsize)/2

  pi=acos(-1.)
  sigma2=2.*sigma**2
  alpha=1/(pi*sigma2)
  
  reduce = 0
  if(wsize.gt.3) then
      reduce = int((wsize-3)/2)
!      mi=mod(reduce,2)
!      print *,reduce
  endif

  do j =1,jsize2
     do i =1,isize2
        sum=0
        mx=xrate*(real(i)-0.5)
        my=yrate*(real(j)-0.5)
        do jj =1+reduce,wsize-reduce
           do ii=1+reduce,wsize-reduce
              i2=int(mx-arad)+ii
              j2=int(my-arad)+jj
              dx=(real(i2)-0.5-mx)
              dy=(real(j2)-0.5-my)
              kernel=alpha*exp(-((dx)**2+(dy)**2)/(sigma2))
              if((i2.le.0).or.(i2.gt.isize)) cycle
              if((j2.le.0).or.(j2.gt.jsize)) cycle
              do k =1,band
                 output(i,j,k)=input(i2,j2,k)*kernel+output(i,j,k)
              enddo
              sum=kernel+sum      
           enddo
        enddo
        do k =1,band
        enddo
        do k =1,band
           output(i,j,k)=output(i,j,k)/sum
        enddo

     enddo
  enddo


END SUBROUTINE gaussian_resize


END MODULE similarity
