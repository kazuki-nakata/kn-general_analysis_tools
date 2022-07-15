MODULE Filter
  IMPLICIT NONE
  INTEGER,PARAMETER :: null = 999
  REAL(4),PARAMETER :: Undef = 999
CONTAINS
  
  SUBROUTINE Index(input,output,isize,jsize,band,band1,band2,mask)
    IMPLICIT NONE
  INTEGER :: i,j,k
  INTEGER(4),INTENT(IN) :: jsize,isize,band,band1,band2,mask
  REAL(4),INTENT(OUT):: output(1:isize,1:jsize)
  REAL(4),INTENT(IN) :: input(1:isize,1:jsize,band)

  do j =1,jsize
     do i =1,isize
!        print *,i,j,input(i,j,band2),input(i,j,band1)
        if(input(i,j,1).eq.mask) then
           output(i,j)=999; cycle
        endif
        
        output(i,j)=(input(i,j,band2)-input(i,j,band1))/(input(i,j,band2)+input(i,j,band1))
        if(output(i,j).ge.1.) output(i,j)=999
        if(output(i,j).le.-1.) output(i,j)=999

     enddo
  enddo
END SUBROUTINE Index

  SUBROUTINE MASK_MERGE(input1,input2,output,isize,jsize)
    IMPLICIT NONE
  INTEGER :: i,j,k
  INTEGER(4),INTENT(IN) :: jsize,isize
  INTEGER(4),INTENT(OUT):: output(1:isize,1:jsize)
  INTEGER(4),INTENT(IN) :: input2(1:isize,1:jsize)
  INTEGER(4),INTENT(IN) ::input1(1:isize,1:jsize)
  output(:,:)=0.
  
  do j =1,jsize
     do i =1,isize
        if((input1(i,j).eq.1).and.(input2(i,j).eq.1)) then
           output(i,j)=1
        else
           output(i,j)=999
        endif
     enddo
  enddo
  
END SUBROUTINE MASK_MERGE


  SUBROUTINE MAKE_MASK_FORLANDSLIDE(input,output,isize,jsize)
    IMPLICIT NONE
  INTEGER :: i,j,k
  INTEGER(4),INTENT(IN) :: jsize,isize
  INTEGER(4),INTENT(OUT):: output(1:isize,1:jsize)
  REAL(4),INTENT(IN) :: input(1:isize,1:jsize,2)

  output(:,:)=1.
  
  do j =1,jsize
     do i =1,isize
        if(input(i,j,1).eq.0) output(i,j)=999.
        if((input(i,j,2).lt.5).or.(input(i,j,2).gt.9)) output(i,j)=999.
     enddo
  enddo
  
END SUBROUTINE MAKE_MASK_FORLANDSLIDE


  SUBROUTINE Interpolation_type1(input,isize,jsize,mask)
    IMPLICIT NONE
  INTEGER :: i,j,k,ii,jj,i2,j2,i3,j3,i4,j4
  INTEGER(4),INTENT(IN) :: jsize,isize
  REAL(4),INTENT(IN) :: mask
  REAL(4),INTENT(INOUT) :: input(1:isize,1:jsize)
  
  do j =1,jsize
     do i =1,isize
! 
        if(input(i,j).eq.mask) then

           if((i-1.lt.1).or.(i+1.gt.isize).or.(j-1.lt.1).or.(j+1.gt.jsize)) then
              if(input(i,j).eq.999) then
                 cycle
              endif
           endif
           
           i3=0;i4=0
           do ii=1,3
              i2=i+ii-4
              if(input(i2,j).eq.mask) then
                 i3=i2
              endif
           enddo

           do ii=1,3
              i2=i+4-ii
              if(input(i2,j).eq.mask) then
                 i4=i2
              endif
           enddo

           if((i3.ne.0).and.(i4.ne.0)) then 
              input(i,j)=(input(i3,j)+input(i4,j))/2.              
           endif
              
           i3=0;i4=0
           do jj=1,3
              j2=j+jj-4
              if(input(i,j2).eq.mask) then
                 j3=j2
              endif
           enddo

           do jj=1,3
              j2=j+4-jj
              if(input(i,j2).eq.mask) then
                 j4=j2
              endif
           enddo

           if((j3.ne.0).and.(j4.ne.0)) then 
              input(i,j)=(input(i,j3)+input(i,j4))/2.              
           endif
           
        endif

     enddo
  enddo
  
END SUBROUTINE Interpolation_Type1

  SUBROUTINE MAKE_MASK_NODATA(input,output,isize,jsize,mask)
    IMPLICIT NONE
  INTEGER :: i,j,k
  INTEGER(4),INTENT(IN) :: jsize,isize
  INTEGER(4),INTENT(IN) :: mask
  INTEGER(4),INTENT(OUT):: output(1:isize,1:jsize)
  REAL(4),INTENT(IN) :: input(1:isize,1:jsize)

  output(:,:)=1.
  
  do j =1,jsize
     do i =1,isize
        if(input(i,j).eq.mask) output(i,j)=999.
     enddo
  enddo
  
END SUBROUTINE MAKE_MASK_NODATA

SUBROUTINE Morphology_Binary_Dilation(input,output,isize,jsize,wsize,mask)
  IMPLICIT NONE
  INTEGER :: i,j,i2,j2,i3,j3,ii,jj,k,mi,mj
  INTEGER(4),INTENT(IN) :: jsize,isize,wsize
  INTEGER(4),INTENT(OUT):: output(1:isize,1:jsize)
  INTEGER(4),INTENT(IN) :: input(1:isize,1:jsize)
  INTEGER(4),ALLOCATABLE :: corse(:,:)
  INTEGER(4) :: isize2,jsize2,mask,n_index
  REAL(4) :: dmi,dmj
  INTEGER(4),PARAMETER :: fore_number=1
  INTEGER(4),PARAMETER :: back_number=0

  mi=mod(isize,wsize)
  mj=mod(jsize,wsize)
  isize2=(isize-mi)/wsize
  jsize2=(jsize-mj)/wsize
  k=(wsize-1)/2

  allocate(corse(-1:isize2+2,-1:jsize2+2))

  corse(-1:isize2+2,-1:jsize2+2)=back_number


  do j =1,jsize2
     do i =1,isize2

        i2=i*wsize-k
        j2=j*wsize-k
        if(input(i2,j2).eq.mask) then
           corse(i,j)=null; cycle
        elseif(input(i2,j2).eq.fore_number) then
           corse(i,j)=fore_number; cycle
        endif

        do jj=1,wsize
           do ii=1,wsize
              i3=i2+ii-k-1
              j3=j2+jj-k-1

              if(input(i3,j3).eq.fore_number) then
                 corse(i,j)=fore_number; exit
              endif
           enddo

           if(input(i3,j3).eq.fore_number) exit
        enddo

     enddo
  enddo

  output(:,:)=0
!  print *,'hi'
  do j =1,jsize
     do i =1,isize
        if(input(i,j).eq.mask) then
           output(i,j)=null; cycle
        elseif(input(i,j).eq.fore_number) then
           output(i,j)=fore_number; cycle
        endif

        mi=mod(i+k,wsize)
        mj=mod(j+k,wsize)
        i2=(i+k-mi)/wsize
        j2=(j+k-mj)/wsize         
        dmi=real(i+k)/real(wsize)
        dmj=mod(j+k,wsize)

        if(dmi-real(mi).ge.0) then
           i3=1
        else
           i3=-1
        endif
        if(dmj-real(mj).ge.0) then
           j3=1
        else
           j3=-1
        endif

        n_index=corse(i2,j2)+corse(i2+i3,j2+j3)+corse(i2,j2+j3)+corse(i2+i3,j2)
        if(n_index.eq.back_number*4) then
           output(i,j)=back_number; cycle
        endif

        output(i,j)=back_number
!        print *,i,j,i3,j3
        do jj=1,wsize
           do ii=1,wsize
              i3=i+ii-k-1
              j3=j+jj-k-1
              if((i3.le.0).or.(i3.gt.isize)) cycle
              if((j3.le.0).or.(j3.gt.jsize)) cycle              
              if(input(i3,j3).eq.fore_number) then
                 output(i,j)=fore_number; exit
              endif
           enddo
        enddo
!        print *,i,j
     enddo
  enddo

END SUBROUTINE Morphology_Binary_Dilation

SUBROUTINE Morphology_Binary_Erosion(input,output,isize,jsize,wsize,mask)
  IMPLICIT NONE
  INTEGER :: i,j,i2,j2,i3,j3,ii,jj,k,mi,mj
  INTEGER(4),INTENT(IN) :: jsize,isize,wsize
  INTEGER(4),INTENT(OUT):: output(1:isize,1:jsize)
  INTEGER(4),INTENT(IN) :: input(1:isize,1:jsize)
  INTEGER(4),ALLOCATABLE :: corse(:,:)
  INTEGER(4) :: isize2,jsize2,mask,n_index
  REAL(4) :: dmi,dmj
  INTEGER(4),PARAMETER :: fore_number=0
  INTEGER(4),PARAMETER :: back_number=1

  mi=mod(isize,wsize)
  mj=mod(jsize,wsize)
  isize2=(isize-mi)/wsize
  jsize2=(jsize-mj)/wsize
  k=(wsize-1)/2

  allocate(corse(-1:isize2+2,-1:jsize2+2))

  corse(-1:isize2+2,-1:jsize2+2)=back_number


  do j =1,jsize2
     do i =1,isize2

        i2=i*wsize-k
        j2=j*wsize-k
        if(input(i2,j2).eq.mask) then
           corse(i,j)=null; cycle
        elseif(input(i2,j2).eq.fore_number) then
           corse(i,j)=fore_number; cycle
        endif

        do jj=1,wsize
           do ii=1,wsize
              i3=i2+ii-k-1
              j3=j2+jj-k-1

              if(input(i3,j3).eq.fore_number) then
                 corse(i,j)=fore_number; exit
              endif
           enddo

           if(input(i3,j3).eq.fore_number) exit
        enddo

     enddo
  enddo

  output(:,:)=0
!  print *,'hi'
  do j =1,jsize
     do i =1,isize    
        if(input(i,j).eq.mask) then
           output(i,j)=null; cycle
        elseif(input(i,j).eq.fore_number) then
           output(i,j)=fore_number; cycle
        endif

        mi=mod(i+k,wsize)
        mj=mod(j+k,wsize)
        i2=(i+k-mi)/wsize
        j2=(j+k-mj)/wsize         
        dmi=real(i+k)/real(wsize)
        dmj=mod(j+k,wsize)

        if(dmi-real(mi).ge.0) then
           i3=1
        else
           i3=-1
        endif
        if(dmj-real(mj).ge.0) then
           j3=1
        else
           j3=-1
        endif
        n_index=corse(i2,j2)+corse(i2+i3,j2+j3)+corse(i2,j2+j3)+corse(i2+i3,j2)
        if(n_index.eq.back_number*4) then
           output(i,j)=back_number; cycle
        endif

        output(i,j)=back_number
!        print *,i,j,i3,j3
        do jj=1,wsize
           do ii=1,wsize
              i3=i+ii-k-1
              j3=j+jj-k-1
              if((i3.le.0).or.(i3.gt.isize)) cycle
              if((j3.le.0).or.(j3.gt.jsize)) cycle              
              if(input(i3,j3).eq.fore_number) then
                 output(i,j)=fore_number; exit
              endif
           enddo
        enddo
!        print *,i,j
     enddo
  enddo

END SUBROUTINE Morphology_Binary_Erosion

  SUBROUTINE Gaussian_filter(input,output,mask,isize,jsize,wsize,sigma,sf)
    IMPLICIT NONE
  INTEGER :: i,j,k,i2,j2,ii,jj
  INTEGER(4),INTENT(IN) :: isize,jsize, wsize
  REAL(4),INTENT(OUT):: output(1:isize,1:jsize)
  REAL(4),INTENT(IN) :: input(isize,jsize)
  INTEGER(4),INTENT(IN) ::mask(isize,jsize)
  REAL(4) :: kernel(wsize,wsize)
  REAL(4) :: sigma,sum,sf

  output(:,:)=0

  k=(wsize-1)/2
  print *,isize,jsize,wsize,sigma
  call gaussian_kernel(kernel,wsize,sigma,sum,sf)

  do j =1,jsize
     do i =1,isize
        if(mask(i,j).eq.999) cycle
        sum=0
        do jj =1,wsize
           do ii=1,wsize
              i2=i+ii-k-1
              j2=j+jj-k-1
              if((i2.le.0).or.(i2.gt.isize)) cycle
              if((j2.le.0).or.(j2.gt.jsize)) cycle

              if(mask(i2,j2).eq.999) cycle

!              print *,output(i,j),input(i2,j2),kernel(ii,jj)
              output(i,j)=input(i2,j2)*kernel(ii,jj)+output(i,j)
              sum=kernel(ii,jj)+sum
        enddo
     enddo
              output(i,j)=output(i,j)/sum
     enddo
  enddo

END SUBROUTINE Gaussian_filter

  SUBROUTINE Binarization(input,output,isize,jsize,threshold,mask)
    IMPLICIT NONE
    INTEGER :: i,j,k,mask
  REAL(4),INTENT(IN) :: threshold
  INTEGER(4),INTENT(IN) :: jsize,isize
  INTEGER(4),INTENT(OUT):: output(1:isize,1:jsize)
  REAL(4),INTENT(IN) :: input(isize,jsize)
  
  do j =1,jsize
     do i =1,isize
        if(input(i,j).eq.mask) then
           output(i,j)=null
        elseif(input(i,j).ge.threshold) then
           output(i,j)=0
        else
           output(i,j)=1
        endif
!        print *,output(i,j),null,mask,threshold,input(i,j)
     enddo
  enddo
END SUBROUTINE Binarization

  SUBROUTINE Raster_Difference_Int(input1,input2,output,isize,jsize,mask)
    IMPLICIT NONE
    INTEGER :: i,j,k,mask
  INTEGER(4),INTENT(IN) :: jsize,isize
  INTEGER(4),INTENT(OUT):: output(1:isize,1:jsize)
  INTEGER(4),INTENT(IN) :: input1(isize,jsize),input2(isize,jsize)
  
  do j =1,jsize
     do i =1,isize
        if(input1(i,j).eq.mask) then
           output(i,j)=mask; cycle
        endif
        
        output(i,j)=input1(i,j)-input2(i,j)
        
     enddo
  enddo
END SUBROUTINE Raster_Difference_Int

  SUBROUTINE down_scaling_real(input,output,isize,jsize,isize2,jsize2,interval)
    IMPLICIT NONE
    INTEGER :: i,j,k,i2,j2,interval
    INTEGER :: isize,jsize,isize2,jsize2
    REAL(4),INTENT(IN) :: input(isize,jsize)
    REAL(4),INTENT(OUT):: output(1:isize2,1:jsize2)

  do j =1,jsize2
     do i =1,isize2
        i2=1+(i-1)*interval
        j2=1+(j-1)*interval
!        if(input(i2,j2).eq.1).or.(input(i2+1,j2).eq.999)+input(i2,j2+1)+input(i2+1,j2+1)
           output(i,j)=(input(i2,j2)+input(i2+1,j2)+input(i2,j2+1)+input(i2+1,j2+1))/4
     enddo
  enddo

END SUBROUTINE down_scaling_real
  SUBROUTINE down_scaling_int(input,output,isize,jsize,isize2,jsize2,interval)
    IMPLICIT NONE
    INTEGER :: i,j,k,i2,j2,interval
    INTEGER :: isize,jsize,isize2,jsize2
    INTEGER(4),INTENT(IN) :: input(isize,jsize)
    INTEGER(4),INTENT(OUT):: output(1:isize2,1:jsize2)

  do j =1,jsize2
     do i =1,isize2
        i2=1+(i-1)*interval
        j2=1+(j-1)*interval
           output(i,j)=(input(i2,j2)+input(i2+1,j2)+input(i2,j2+1)+input(i2+1,j2+1))/4
     enddo
  enddo

END SUBROUTINE down_scaling_int


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

SUBROUTINE moving_average(width,length,band,mask,mask2,wsize,data)
   INTEGER,INTENT(IN)::width
   INTEGER,INTENT(IN)::length
   INTEGER,INTENT(IN)::band
   INTEGER,INTENT(IN)::wsize
   REAL(4)::dum(width,length,band)
   REAL(4),INTENT(INOUT) :: data(width,length,band)
   REAL(4),INTENT(IN) :: mask,mask2
   INTEGER::i,j,k,ii,jj,i3,j3,inum
   dum(:,:,:)=0.
   do j =1,length
      do i =1,width
!          if((data(i,j,1).ge.-0.1).and.(data(i,j,1).le.0.1)) then
!          print *,data(i,j,1)
!          endif
         if(data(i,j,1).ne.int(mask)) then
            inum=0
         do jj=1,wsize
            do ii=1,wsize
               i3=i+ii-(wsize-1)/2-1
               j3=j+jj-(wsize-1)/2-1
               if((i3.lt.1).or.(i3.gt.width)) cycle
               if((j3.lt.1).or.(j3.gt.length)) cycle
               if(data(i3,j3,1).eq.mask) cycle
               inum=1+inum            
               do k=1,band            
                  dum(i,j,k)=data(i3,j3,k)+dum(i,j,k)
               enddo
            enddo
         enddo
         dum(i,j,:)=dum(i,j,:)/real(inum)
!          dum(i,j,2)=dum(i,j,2)/real(inum)
!          dum(i,j,3)=dum(i,j,3)/real(inum)
      else
         dum(i,j,:)=mask2
         endif
      enddo
   enddo
   data(:,:,:)=dum(:,:,:)
END SUBROUTINE moving_average

SUBROUTINE Interpolation_type2(inval,isize2,jsize2,inlat,inlon,isize,jsize,ymin,xmin,yint,xint,output)
   IMPLICIT NONE
   INTEGER :: i,j,ii,jj,k
   iNTEGER :: ipos,jpos,ipos2,jpos2
   REAL(4) :: dx,dy
   INTEGER(4),INTENT(IN) :: jsize,isize,isize2,jsize2
   REAL(4),INTENT(IN) :: xmin,ymin,xint,yint
   REAL(4),INTENT(IN) :: inlat(1:isize,1:jsize),inlon(1:isize,1:jsize)
   REAL(4),INTENT(IN) :: inval(1:isize2,1:jsize2)
   REAL(4) :: rlon(1:jsize2),rlat(1:isize2)
   INTEGER(4),INTENT(OUT):: output(1:isize,1:jsize)

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
               ipos=jj
               ipos2=jj-1
            endif
         enddo

         ipos=-999
         do ii = 2, isize2
            if((inlat(i,j) > rlat(ii-1)) .and. &
                  (inlat(i,j) <= rlat(ii))) then
            ipos=ii
            ipos2=ii-1
            endif
         enddo

         if(jpos == -999) then
            jpos2=1
            jpos=jsize2
         endif

         if(jpos == -999) then
            cycle
         endif

         if(ipos == -999) then
         dx=(inlon(i,j)-rlon(jpos2))  / ((rlon(jpos)+360.)-rlon(jpos2))
         dy=(inlat(i,j)-rlat(ipos2))  / (rlat(ipos)-rlat(ipos2))
            else
         dx=(inlon(i,j)-rlon(jpos2)) / (rlon(jpos)-rlon(jpos2))
         dy=(inlat(i,j)-rlat(ipos2))  / (rlat(ipos)-rlat(ipos2))
         endif
         print *,ipos,ipos2,jpos,jpos2,isize2,jsize2
         ! do k = 1, 1!Nvars
         !    output(i,j) = bilin(dy,dx,Undef, &
         !                   inval(ipos2,jpos2), &
         !                   inval(ipos  ,jpos2), &
         !                   inval(ipos  ,jpos ), &
         !                   inval(ipos2,jpos  ))
         ! enddo
      enddo
      
      enddo
      output(:,:)=0.
END SUBROUTINE Interpolation_Type2

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

END MODULE Filter
