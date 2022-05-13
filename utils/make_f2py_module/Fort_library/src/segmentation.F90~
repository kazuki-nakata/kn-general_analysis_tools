MODULE segmentation
  use OBJECT_MODULE
  IMPLICIT NONE
  INTEGER(4) :: isize,jsize,band  
  integer(4) :: i,j,k,kk,i1,j1,i2,j2,i3,j3
  integer(4) :: irec,ii,jj,ui    
  integer(4), allocatable :: pairs1(:),pairs2(:)
  integer(4),parameter :: band_num=3
  REAL(4):: diff(band_num)
  REAL(4), ALLOCATABLE:: segment(:,:,:),data(:,:,:) 
  INTEGER(4),ALLOCATABLE::edge(:,:,:)
  real(4) :: mask,logdelta
  real(4) :: smallregion,red,green,blue,size
  real(4), allocatable::diffmax(:)
  integer(4) :: pmax_uf
  integer(4) :: ipos_8c(4),jpos_8c(4)
  integer(4) :: ipos_4c(2),jpos_4c(2)
  REAL(4) :: q  
  REAL(4),PARAMETER ::thres_size_min=60   ! 36    !-----i.e. 18m x 18m, 12m x 27m
  REAL(4),PARAMETER ::thres_size_max=25000            !---- i.e. 1000m x 250m 
!   REAL(4),PARAMETER ::thres_ls1= 1.5
!   REAL(4),PARAMETER ::thres_ls2= 0.8
!   REAL(4),PARAMETER ::thres_red= 100.
  REAL(4),PARAMETER ::thres_ls1= 1.3
  REAL(4),PARAMETER ::thres_ls2= 0.8
  REAL(4),PARAMETER ::thres_red= 160.  
  REAL(4),PARAMETER ::thres_asp = 0.8 !0.8
  REAL(4),PARAMETER ::thres_comp = 0.35 !0.30 !0.35
  INTEGER(4) :: irec1,irec2
  INTEGER(4) :: dummy1(100000),dummy2(100000)
  INTEGER(4),ALLOCATABLE :: red_object(:),green_object(:)
  REAL(4),ALLOCATABLE :: osize(:),compactness(:),aspectratio(:),symmetry(:),concomplex(:)
  REAL(4) :: imin(13),imax(13),jmin(13),jmax(13)
  REAL(4) :: xmin(13),xmax(13),ymin(13),ymax(13),amin,asp
  REAL(4) :: cos_as(13),sin_as(13),rot_x,rot_y,pi,dora
  INTEGER(4),PARAMETER::null=0
  
CONTAINS
  
  subroutine landslide_detection(data_i,data_o,band2,length,width)
  use OBJECT_MODULE      
    implicit none  
   INTEGER(4),INTENT(IN) :: band2,length,width
   REAL(4),INTENT(IN) :: data_i(band2,length,width)  
   INTEGER(4),INTENT(OUT) :: data_o(length,width)     
  pi=acos(-1.)
  dora=pi/180.
  !-----------------data reading---------------------------------
  isize=width
  jsize=length
  band=band2
  allocate(data(isize,jsize,band))
  allocate(segment(isize,jsize,band))
  allocate(edge(isize,jsize,2))
  print *,isize,jsize,band
  CALL reshape_backforward_3D(data_i(1:band,1:jsize,1:isize),data(1:isize,1:jsize,1:band),band,jsize,isize)
  print *,'hi'  
  call moving_average(isize,jsize,band,0,999,3,data(:,:,:))
  print *,'hi'
  !---------------linear transformation to RGB-------------------
  do j=1,jsize
     do i=1,isize
        if(int(data(i,j,1)).eq.999) cycle
!   data(i,j,1)=255./20.*(25.+data(i,j,1))
!   data(i,j,2)=255./20.*(25.+data(i,j,2))
!   data(i,j,3)=255./10*data(i,j,3)        
   data(i,j,1)=255./20.*(25.+data(i,j,1))
   data(i,j,2)=255./20.*(25.+data(i,j,2))
   data(i,j,3)=255./10*data(i,j,3)
   if(data(i,j,1).gt.255) data(i,j,1)=255
   if(data(i,j,2).gt.255) data(i,j,2)=255
   if(data(i,j,3).gt.255) data(i,j,3)=255
   if(data(i,j,1).lt.1) data(i,j,1)=1
   if(data(i,j,2).lt.1) data(i,j,2)=1
   if(data(i,j,3).lt.1) data(i,j,3)=1   
enddo
enddo

  !---------first segmentation--------------------
q=256
!mode=1
  call statistical_region_merging
  call deallocate_uf

  !-------second segmentation---------------------
  q=128
!  mode=2
  data(:,:,:)=segment(:,:,:)
  call statistical_region_merging
  print *, "finish srm"
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
  do while (ii.ne.smax_uf+1)

     irec=1+irec
     size=get_size_uf(ii)
     red=get_property_uf(ii,1)
     green=get_property_uf(ii,2)
     blue=get_property_uf(ii,3)
     if((size.ge.thres_size_max).or.(red.le.thres_red).or.(red/green.lt.thres_ls1)) then
        ii=get_par_dll(ii);cycle      
     endif
     
     jj=get_nchi_sll(ii)
     print *,ii,jj,ii-jj
     if(red/green.ge.thres_ls1) then
        
        do while (jj.ne.NULL)
           i2=get_ipos_uf(jj);j2=get_jpos_uf(jj)
           i3=edge(i2,j2,1);j3=edge(i2,j2,2)

           if((i3.ne.0).and.(j3.ne.0)) then
              kk=find_uf(ptrt(i3,j3))
              red=get_property_uf(kk,1)
              green=get_property_uf(kk,2)
              if(red/green.ge.thres_ls1) then
                 call union_uf(ii,kk,ui)
              endif
           endif
           jj=get_nchi_sll(jj)
           if(size.eq.1) jj=0
        enddo
        irec1=irec1+1
        dummy1(irec1)=find_uf(ii)
     endif
     
     ii=get_par_dll(ii)
  enddo
  print *,'demo'
!------------------list arrangement and filtering small region------
  irec=0
do i=1,irec1
   ii=dummy1(i)
     size=get_size_uf(ii)
     red=get_property_uf(ii,1)
     if((size.gt.thres_size_min).and.(size.lt.thres_size_max).and.(red.ge.thres_red)) then   !also remove size 0
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

   do while (jj.ne.NULL)
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
!!$   red=get_property_uf(ii,1)
!!$   green=get_property_uf(ii,2)
!!$   blue=get_property_uf(ii,3)
   i2=ii
   if((compactness(i).le.thres_comp).and.(aspectratio(i).le.thres_asp)) then
   do while(i2.ne.NULL)
      data(get_ipos_uf(i2),get_jpos_uf(i2),1)=1
      i2=get_nchi_sll(i2)
   enddo
   endif
   enddo
   
  CALL reshape_backforward_2D(data(:,:,1),data_o(:,:),isize,jsize)     
  !----------------------------------------------
endsubroutine landslide_detection
   
SUBROUTINE statistical_region_merging
  IMPLICIT NONE
   mask=999
  !-------------initialization of uf----------------------------     
   call init_unionfind(isize,jsize,data(:,:,:),band_num,mask,smax_uf)     
  !------------connected system-------------------------------
  ipos_8c=(/-1,0,1,1/)
  jpos_8c=(/1,1,1,0/)
  ipos_4c=(/1,0/)
  jpos_4c=(/0,1/)
  allocate(pairs1(smax_uf*4),pairs2(smax_uf*4))
  pairs1(:)=-999;pairs2(:)=-999
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
!-------------heap sort by gradients between pixel values-------
  allocate(diffmax(1:pmax_uf))
  do i=1,pmax_uf
    do k=1,3
     diff(k)=abs(get_property_uf(pairs1(i),k)-get_property_uf(pairs2(i),k))
!     print *,diff(k),k
  enddo
!  diff(1)=diff(1)/diff(3)
!  diff(2)=diff(2)/diff(3)
  diffmax(i)=maxval(diff(1:3))
!  print *,diffmax(i)
     enddo
!  print *,
     call heapsort(pmax_uf,diffmax,pairs1,pairs2)
!---------------region merging and edge detection---------------------
     irec=0
     edge(:,:,:)=0
     do i=1,pmax_uf
         i2=get_ipos_uf(pairs1(i))
         j2=get_jpos_uf(pairs1(i))
         i3=get_ipos_uf(pairs2(i))
         j3=get_jpos_uf(pairs2(i))        
         !     print *,pairs1(i),pairs2(i)
!         print *,predicate(pairs1(i),pairs2(i)),i
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
   segment(:,:,:)=0.
  do j=1,jsize
     do i=1,isize
        ii=find_uf(ptrt(i,j))
        do k=1,3
           segment(i,j,k)=get_property_uf(ii,k)
        enddo
  enddo
  enddo

  DEALLOCATE(diffmax,pairs1,pairs2)
end subroutine STATISTICAL_REGION_MERGING


LOGICAL FUNCTION predicate(tmp1,tmp2)
  use OBJECT_MODULE
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

    do j =1,3
       diff(j)=(get_property_uf(r1,j)-get_property_uf(r2,j))**2
    enddo
    r1_sz=real(get_size_uf(r1))
    r2_sz=real(get_size_uf(r2))
    logr1=min(g,r1_sz)*log(1.0+r1_sz)
    logr2=min(g,r2_sz)*log(1.0+r2_sz)
    b=g**2./(2.*q*r1_sz)*(logr1+logdelta)
    b=b+g**2./(2.*q*r2_sz)*(logr2+logdelta)
!    if((diff(1).le.b).and.(diff(2).le.b).and.(diff(3).le.b)) then
    if((diff(1).le.b).and.(diff(2).le.b)) then    
       predicate=.true.
       return
    else
       predicate=.false.
    endif
  END FUNCTION predicate

END MODULE Segmentation
