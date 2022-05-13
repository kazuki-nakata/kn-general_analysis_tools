MODULE segmentation_RB
  use OBJECT_MODULE2
  IMPLICIT NONE
  INTEGER(4) :: isize,jsize,band  
  integer(4) :: i,j,k,kk,i1,j1,i2,j2,i3,j3
  integer(4) :: irec,ii,jj,ui    
  integer(4), allocatable :: pairs1(:),pairs2(:)
  integer(4),parameter :: band_num=2
  REAL(4):: diff(band_num)
  REAL(4), ALLOCATABLE:: segment(:,:,:),data(:,:,:) 
  INTEGER(4),ALLOCATABLE::edge(:,:,:)
  real(4) :: mask,logdelta
  real(4) :: smallregion,red,green,blue,size,black_object
  real(4), allocatable::diffmax(:)
  integer(4) :: pmax_uf
  integer(4) :: ipos_8c(4),jpos_8c(4)
  integer(4) :: ipos_4c(2),jpos_4c(2)
  REAL(4) :: q  
  INTEGER(4) :: irec1,irec2
  INTEGER(4) :: dummy1(1000000),dummy2(1000000)
  INTEGER(4),ALLOCATABLE :: red_object(:),blue_object(:)
  REAL(4),ALLOCATABLE :: osize(:),compactness(:),aspectratio(:),symmetry(:),concomplex(:)
  REAL(4) :: imin(13),imax(13),jmin(13),jmax(13)
  REAL(4) :: xmin(13),xmax(13),ymin(13),ymax(13),amin,asp
  REAL(4) :: cos_as(13),sin_as(13),rot_x,rot_y,pi,dora
  INTEGER(4),PARAMETER::iNULL=0
  
CONTAINS
  
  subroutine landslide_detection(data_i,data_o,band2,length,width,sigma_thres1,sigma_thres2, &
       color_thres,comp_thres,asp_thres,size_min_thres)
    use OBJECT_MODULE2
    use Basic_Processing
    implicit none  
   INTEGER(4),INTENT(IN) :: band2,length,width
   REAL(4),INTENT(IN) :: data_i(band2,length,width)  
   INTEGER(4),INTENT(OUT) :: data_o(length,width)
   REAL(4),INTENT(IN) :: sigma_thres1
   REAL(4),INTENT(IN) :: sigma_thres2
   REAL(4),INTENT(IN) :: color_thres
   REAL(4),INTENT(IN) :: comp_thres
   REAL(4),INTENT(IN) :: asp_thres
   REAL(4),INTENT(IN) :: size_min_thres
  REAL(4),PARAMETER ::thres_size_max=30000            !---- i.e. 1000m x 300m 
  REAL(4),PARAMETER ::thres_red= 80.  
  INTEGER(4),ALLOCATABLE :: mask(:,:),mask2(:,:)
  !------test data---
  REAL(4), ALLOCATABLE :: test(:,:,:)
  character(512) :: dir
  !---------  
  pi=acos(-1.)
  dora=pi/180.
  !-----------------data reading---------------------------------
  isize=width
  jsize=length
  band=band2
  allocate(test(isize,jsize,3))   
  allocate(data(isize,jsize,band))
  allocate(segment(isize,jsize,band))
  allocate(edge(isize,jsize,2))
  allocate(mask(isize,jsize),mask2(isize,jsize))
  print *,isize,jsize,band
  CALL reshape_backforward_3D(data_i(1:band,1:jsize,1:isize),data(1:isize,1:jsize,1:band),band,jsize,isize)

  call moving_average(isize,jsize,band,0,999,3,data(:,:,:))
 

  mask(:,:)=0
  !----------------------data mask--------------------------------
  do j=1,jsize
     do i=1,isize
        if(data(i,j,1).eq.0) then
           data(i,j,:)=999;mask(i,j)=999
        elseif(data(i,j,2).eq.0) then
           data(i,j,:)=999;mask(i,j)=999
        endif
        if(data(i,j,1).ge.sigma_thres1) mask(i,j)=1
        if(data(i,j,2).ge.sigma_thres1) mask(i,j)=1
     enddo
  enddo
  
  !--------------MAKE MASK DATA FOR SMALL SCALE REGION----------
  DO i =1,1
     CALL Morphology_Binary_Erosion(mask(:,:),mask2(:,:),isize,jsize,3,999)
     mask(:,:)=mask2(:,:)
  ENDDO  
  DO i =1,4
     CALL Morphology_Binary_Dilation(mask(:,:),mask2(:,:),isize,jsize,3,999)
     mask(:,:)=mask2(:,:)
  ENDDO

  do j=1,jsize
     do i=1,isize
        if(mask(i,j).eq.1) data(i,j,:)=999
     enddo
  enddo
  
  deallocate(mask,mask2)
  
  !---------------linear transformation to RGB-------------------
  do j=1,jsize
     do i=1,isize
        if(int(data(i,j,1)).eq.999) cycle   
        data(i,j,1)=255./20.*(data(i,j,1)-60)
        data(i,j,2)=255./20.*(data(i,j,2)-60)
        if(data(i,j,1).gt.255) data(i,j,1)=255
        if(data(i,j,2).gt.255) data(i,j,2)=255
        if(data(i,j,1).lt.1) data(i,j,1)=1
        if(data(i,j,2).lt.1) data(i,j,2)=1
     enddo
  enddo

  black_object=255./20. * (sigma_thres2-60)
!--------------------test---------------------------------
!!$  dir='D:\segmentation\FY30\\mask.tif'
!!$  test(:,:,1)=data(:,:,1)
!!$  test(:,:,2)=data(:,:,2)
!!$  test(:,:,3)=data(:,:,2)
!!$  call write_tif(dir,test,isize,jsize,3,jsize)
  
  !---------first segmentation--------------------
q=256
  call statistical_region_merging
  call deallocate_uf

  !-------second segmentation---------------------
q=256
  data(:,:,:)=segment(:,:,:)
  call statistical_region_merging
!  print *, "finish srm"
  
!--------------------test---------------------------------
!!$  dir='D:\segmentation\FY30\\segment3.tif'
!!$  test(:,:,1)=segment(:,:,1)
!!$  test(:,:,2)=segment(:,:,2)
!!$  test(:,:,3)=segment(:,:,2)
!!$  call write_tif(dir,test,isize,jsize,3,jsize)
  
  !---------identification of first child node --------
  irec=0
  do i =1, smax_uf
     if(get_chi_dll(i).eq.0) then
        ii=i;exit
     endif
  enddo
  irec1=0
  irec2=0

  !----------merging between residual similar color regions-------
  print *, '------merging between residual similar color regions-----'
  do while (ii.ne.smax_uf+1)
     
     irec=1+irec
     size=get_size_uf(ii)
     red=get_property_uf(ii,1)
     blue=get_property_uf(ii,2)
     if((size.ge.thres_size_max).or.(red.le.thres_red).or.(red/blue.lt.color_thres).or.((red+blue)/2.lt.black_object)) then
        ii=get_par_dll(ii);cycle      
     endif
     
     jj=get_nchi_sll(ii)
     if(red/blue.ge.color_thres) then
        
        do while (jj.ne.iNULL)
           i2=get_ipos_uf(jj);j2=get_jpos_uf(jj)
           i3=edge(i2,j2,1);j3=edge(i2,j2,2)
           if((i3.ne.0).and.(j3.ne.0)) then
              kk=find_uf(ptrt(i3,j3))
              red=get_property_uf(kk,1)
              blue=get_property_uf(kk,2)        
              if((red/blue.ge.color_thres).and.((red+blue)/2.ge.black_object)) then
                 call union_uf(ii,kk,ui)
              endif
           endif
           jj=get_nchi_sll(jj)
           if(size.eq.1) jj=0
!           print *,jj,ii
        enddo
        irec1=irec1+1
        dummy1(irec1)=find_uf(ii)
     endif
     
     ii=get_par_dll(ii)
  enddo

  !------------------list arrangement and filtering small region------
  print *, '---list arrangement and filtering small region---'
  irec=0
do i=1,irec1
   ii=dummy1(i)
     size=get_size_uf(ii)
     red=get_property_uf(ii,1)
     if((size.gt.size_min_thres).and.(size.lt.thres_size_max).and.(red.ge.thres_red)) then   !also remove size 0
        irec=1+irec
     dummy2(irec)=ii
  endif
enddo
allocate(red_object(1:irec))
red_object(1:irec)=dummy2(1:irec)
print *,'size and shape'
!--------estimation of natural shape statistics--------
allocate(osize(1:irec),compactness(1:irec),aspectratio(1:irec),symmetry(1:irec),concomplex(1:irec))
compactness(:)=0.
do i=1,13
   cos_as(i)=cos(real(i-1)*15*dora)
   sin_as(i)=sin(real(i-1)*15*dora)
enddo

do i=1,irec
   
   ii=red_object(i)
   imin(:)=1000000;imax(:)=-1000000;jmin(:)=1000000;jmax(:)=-1000000;   
   osize(i)=get_size_uf(ii)
   jj=get_nchi_sll(ii)

   do while (jj.ne.iNULL)
      i2=get_ipos_uf(jj)
      j2=get_jpos_uf(jj)
!--------------compactness---------------------------
      if(edge(i2,j2,1).ne.0) then
         compactness(i)=compactness(i)+1.
      endif
!------------------aspect ratio------------------
      do k = 1, 13
         rot_x=real(i2)*cos_as(k)+real(j2)*sin_as(k)
         rot_y=real(j2)*cos_as(k)-real(i2)*sin_as(k)
         imin(k)=min(rot_x,imin(k))
         if(imin(k).eq.rot_x) xmin(k)=real(imin(k))
         imax(k)=max(rot_x,imax(k))
         if(imax(k).eq.rot_x) xmax(k)=real(imax(k))
         jmin(k)=min(rot_y,jmin(k))
         if(jmin(k).eq.rot_y) ymin(k)=real(jmin(k))
         jmax(k)=max(rot_y,jmax(k))
         if(jmax(k).eq.rot_y) ymax(k)=real(jmax(k))
      enddo
 
      jj=get_nchi_sll(jj)
   enddo

   amin=1000000        
   do k =1,13
      asp=abs(xmax(k)-xmin(k))/abs(ymax(k)-ymin(k))
      amin=min(asp,amin)
      if(amin.eq.asp) aspectratio(i)=amin
   enddo
   compactness(i)=compactness(i)/osize(i)
enddo



data(:,:,:)=0
do i =1, irec
   ii=red_object(i)
!   red=get_property_uf(ii,1)
!   blue=get_property_uf(ii,2)
   !   blue=get_property_uf(ii,3)
!   print *, compactness(i),aspectratio(i)
   i2=ii
   if((compactness(i).le.comp_thres).and.(aspectratio(i).le.asp_thres)) then
   do while(i2.ne.iNULL)
      data(get_ipos_uf(i2),get_jpos_uf(i2),1)=1
      i2=get_nchi_sll(i2)
   enddo
   endif
   enddo

!!$  dir='D:\segmentation\FY30\\result.tif'
!!$  test(:,:,1)=segment(:,:,1)
!!$  test(:,:,2)=segment(:,:,2)
!!$  test(:,:,3)=data(:,:,1)
!!$  call write_tif(dir,test,isize,jsize,3,jsize)
   
  CALL reshape_backforward_2D(data(:,:,1),data_o(:,:),isize,jsize)     
  !----------------------------------------------
endsubroutine landslide_detection
   
SUBROUTINE statistical_region_merging
  IMPLICIT NONE
  mask=999
  print *,'-----------statistical_region_merging---------------'
  !-------------initialization of uf----------------------------
  print *, 'initialize!'
  call init_unionfind(isize,jsize,data(:,:,:),band_num,mask,smax_uf)     
  !------------connected system-------------------------------
  ipos_8c=(/-1,0,1,1/)
  jpos_8c=(/1,1,1,0/)
  ipos_4c=(/1,0/)
  jpos_4c=(/0,1/)
  allocate(pairs1(smax_uf*4),pairs2(smax_uf*4))
  pairs1(:)=-999;pairs2(:)=-999

  print *, 'make pixel pairs'
  irec=0
  do j=1,jsize
     do i=1,isize
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
  !-------------------------------------------------------------------  
  pmax_uf=irec
  logdelta = 2.0 * log(6.0 * smax_uf)
  logdelta = 2.0 * log(6.0 * smax_uf/10000)
  !-------------heap sort by gradients between pixel values-------
  allocate(diffmax(1:pmax_uf))
  do i=1,pmax_uf
!     do k=1,2
!        diff(k)=abs(get_property_uf(pairs1(i),k)-get_property_uf(pairs2(i),k))
!     enddo
     diff(1)=get_property_uf(pairs1(i),1)-get_property_uf(pairs1(i),2)
     diff(2)=get_property_uf(pairs2(i),1)-get_property_uf(pairs2(i),2)
     diffmax(i)=abs(diff(1)-diff(2))
!     diffmax(i)=maxval(diff(1:2))
  enddo

  print *,'sort for horizontal differences'
  call heapsort(pmax_uf,diffmax,pairs1,pairs2)
  !---------------region merging and edge detection---------------------

  print *,'union based on a predicate. pmax_uf=', pmax_uf/1000/1000, 'M'
  irec=0
  edge(:,:,:)=0
  do i=1,pmax_uf
     i2=get_ipos_uf(pairs1(i))
     j2=get_jpos_uf(pairs1(i))
     i3=get_ipos_uf(pairs2(i))
     j3=get_jpos_uf(pairs2(i))        

     if(predicate(pairs1(i),pairs2(i))) then
        call union_uf(pairs1(i),pairs2(i),ui)
        edge(i2,j2,1)=0
        edge(i2,j2,2)=0
     else
        if(edge(i3,j3,1).eq.0) then
           edge(i2,j2,1)=i3
           edge(i2,j2,2)=j3
        endif
     endif
  enddo
  print *, 'make a segment image'
  segment(:,:,:)=999.
  do j=1,jsize
     do i=1,isize
        if(data(i,j,1).eq.999) cycle
!        print *,i,j,data(i,j,1),data(i,j,2)
        ii=find_uf(ptrt(i,j))
!        print *,i,j,ii,data(i,j,1),data(i,j,2)
        do k=1,2
           segment(i,j,k)=get_property_uf(ii,k)
        enddo
     enddo
  enddo

  DEALLOCATE(diffmax,pairs1,pairs2)
end subroutine STATISTICAL_REGION_MERGING


LOGICAL FUNCTION predicate(tmp1,tmp2)
  use OBJECT_MODULE2
  implicit none
  INTEGER(4),INTENT(IN)::tmp1,tmp2
  INTEGER(4)::r1,r2
  REAL(4),PARAMETER :: g =256.0
  REAL(4) :: b
  REAL(4) :: r1_sz,r2_sz,logr1,logr2
  REAL(4) :: diff(3)
  !    print *,pmax_uf

  r1=find_uf(tmp1)
  r2=find_uf(tmp2)

    do j =1,2
       diff(j)=(get_property_uf(r1,j)-get_property_uf(r2,j))**2
    enddo
    r1_sz=real(get_size_uf(r1))
    r2_sz=real(get_size_uf(r2))
    logr1=min(g,r1_sz)*log(1.0+r1_sz)
    logr2=min(g,r2_sz)*log(1.0+r2_sz)
    b=g**2./(2.*q*r1_sz)*(logr1+logdelta)
    b=b+g**2./(2.*q*r2_sz)*(logr2+logdelta)
    if((diff(1).le.b).and.(diff(2).le.b)) then    
       predicate=.true.
       return
    else
       predicate=.false.
    endif
  END FUNCTION predicate


SUBROUTINE calculate_average_and_std(data_i,band,jsize,isize,ave,std)
  IMPLICIT NONE
  INTEGER :: i,j,k,ii,jj
  REAL(4):: xrate,yrate,xoffs,yoffs,num
  INTEGER(4),INTENT(IN) :: isize,jsize,band
  REAL(4), INTENT(OUT)  :: ave(2),std(2)
  REAL(4), INTENT(IN) :: data_i(band,jsize,isize)  
  
  ave(:)=0
  std(:)=0
  do k =1,band
     num = 0
     do j =1,jsize
        do i=1,isize
           if(data_i(k,j,i).eq.0) then
              ave(k)=data_i(k,j,i)+ave(k)
              num = num + 1
           endif
        enddo
     enddo
     ave(k)=ave(k)/num
  enddo

  do k =1,band
     num = 0
     do j =1,jsize
        do i=1,isize
           if(data_i(k,j,i).eq.0) then
              std(k)=(data_i(k,j,i)-ave(k))**2+std(k)
              num = num + 1
           endif
        enddo
     enddo
     std(k)=sqrt(std(k)/num)
  enddo
  
END SUBROUTINE calculate_average_and_std
  
END MODULE Segmentation_RB
